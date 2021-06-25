//*********************************************************
//
// Copyright (c) Microsoft. All rights reserved.
// This code is licensed under the MIT License (MIT).
// THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
// IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
// PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
//
//*********************************************************

#ifndef RAYTRACING_HLSL
#define RAYTRACING_HLSL

#define HLSL
#include "RaytracingHlslCompat.h"
#include "ProceduralPrimitivesLibrary.hlsli"
#include "RaytracingShaderHelper.hlsli"
static const float Pi = 3.141592654f;
static const float Pi2 = 6.283185307f;
static const float Pi_2 = 1.570796327f;
static const float Pi_4 = 0.7853981635f;
static const float InvPi = 0.318309886f;
static const float InvPi2 = 0.159154943f;
//***************************************************************************
//*****------ Shader resources bound via root signatures -------*************
//***************************************************************************

// Scene wide resources.
//  g_* - bound via a global root signature.
//  l_* - bound via a local root signature.
RaytracingAccelerationStructure g_scene : register(t0, space0);
RWTexture2D<float4> g_renderTarget : register(u0);
ConstantBuffer<SceneConstantBuffer> g_sceneCB : register(b0);

// Triangle resources
ByteAddressBuffer g_indices : register(t1, space0);
StructuredBuffer<Vertex> g_vertices : register(t2, space0);

// Procedural geometry resources
StructuredBuffer<PrimitiveInstancePerFrameBuffer> g_AABBPrimitiveAttributes : register(t3, space0);
ConstantBuffer<PrimitiveConstantBuffer> l_materialCB : register(b1);
ConstantBuffer<PrimitiveInstanceConstantBuffer> l_aabbCB: register(b2);


//***************************************************************************
//****************------ Utility functions -------***************************
//***************************************************************************

// Diffuse lighting calculation.
float CalculateDiffuseCoefficient(in float3 hitPosition, in float3 incidentLightRay, in float3 normal)
{
    float fNDotL = saturate(dot(-incidentLightRay, normal));
    return fNDotL;
}

// Phong lighting specular component
float4 CalculateSpecularCoefficient(in float3 hitPosition, in float3 incidentLightRay, in float3 normal, in float specularPower)
{
    float3 reflectedLightRay = normalize(reflect(incidentLightRay, normal));
    return pow(saturate(dot(reflectedLightRay, normalize (-WorldRayDirection()))), specularPower);
}


// Phong lighting model = ambient + diffuse + specular components.
float4 CalculatePhongLighting(in float4 albedo, in float3 normal, in bool isInShadow, in float diffuseCoef = 1.0, in float specularCoef = 1.0, in float specularPower = 50)
{
    float3 hitPosition = HitWorldPosition();
    float3 lightPosition = g_sceneCB.lightPosition.xyz;
    float shadowFactor = isInShadow ? InShadowRadiance : 1.0;
    float3 incidentLightRay = normalize(hitPosition - lightPosition);

    // Diffuse component.
    float4 lightDiffuseColor = g_sceneCB.lightDiffuseColor;
    float Kd = CalculateDiffuseCoefficient(hitPosition, incidentLightRay, normal);
    float4 diffuseColor = shadowFactor * diffuseCoef * Kd * lightDiffuseColor * albedo;

    // Specular component.
    float4 specularColor = float4(0, 0, 0, 0);
    if (!isInShadow)
    {
        float4 lightSpecularColor = float4(1, 1, 1, 1);
        float4 Ks = CalculateSpecularCoefficient(hitPosition, incidentLightRay, normal, specularPower);
        specularColor = specularCoef * Ks * lightSpecularColor;
    }

    // Ambient component.
    // Fake AO: Darken faces with normal facing downwards/away from the sky a little bit.
    float4 ambientColor = g_sceneCB.lightAmbientColor;
    float4 ambientColorMin = g_sceneCB.lightAmbientColor - 0.1;
    float4 ambientColorMax = g_sceneCB.lightAmbientColor;
    float a = 1 - saturate(dot(normal, float3(0, -1, 0)));
    ambientColor = albedo * lerp(ambientColorMin, ambientColorMax, a);
    
    return ambientColor + diffuseColor + specularColor;
}

//***************************************************************************
//*****------ TraceRay wrappers for radiance and shadow rays. -------********
//***************************************************************************

// Trace a radiance ray into the scene and returns a shaded color.
float4 TraceRadianceRay(in Ray ray, in UINT currentRayRecursionDepth)
{
    if (currentRayRecursionDepth >= MAX_RAY_RECURSION_DEPTH)
    {
        return float4(0, 0, 0, 0);
    }

    // Set the ray's extents.
    RayDesc rayDesc;
    rayDesc.Origin = ray.origin;
    rayDesc.Direction = ray.direction;
    // Set TMin to a zero value to avoid aliasing artifacts along contact areas.
    // Note: make sure to enable face culling so as to avoid surface face fighting.
    rayDesc.TMin = 0;
    rayDesc.TMax = 10000;
    RayPayload rayPayload = { float4(0, 0, 0, 0), currentRayRecursionDepth + 1 };
    TraceRay(g_scene,
        RAY_FLAG_CULL_BACK_FACING_TRIANGLES,
        TraceRayParameters::InstanceMask,
        TraceRayParameters::HitGroup::Offset[RayType::Radiance],
        TraceRayParameters::HitGroup::GeometryStride,
        TraceRayParameters::MissShader::Offset[RayType::Radiance],
        rayDesc, rayPayload);

    return rayPayload.color;
}

// Trace a shadow ray and return true if it hits any geometry.
bool TraceShadowRayAndReportIfHit(in Ray ray, in UINT currentRayRecursionDepth)
{
    if (currentRayRecursionDepth >= MAX_RAY_RECURSION_DEPTH)
    {
        return false;
    }

    // Set the ray's extents.
    RayDesc rayDesc;
    rayDesc.Origin = ray.origin;
    rayDesc.Direction = ray.direction;
    // Set TMin to a zero value to avoid aliasing artifcats along contact areas.
    // Note: make sure to enable back-face culling so as to avoid surface face fighting.
    rayDesc.TMin = 0;
    rayDesc.TMax = 10000;

    // Initialize shadow ray payload.
    // Set the initial value to true since closest and any hit shaders are skipped. 
    // Shadow miss shader, if called, will set it to false.
    ShadowRayPayload shadowPayload = { true };
    TraceRay(g_scene,
        RAY_FLAG_CULL_BACK_FACING_TRIANGLES
        | RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH
        | RAY_FLAG_FORCE_OPAQUE             // ~skip any hit shaders
        | RAY_FLAG_SKIP_CLOSEST_HIT_SHADER, // ~skip closest hit shaders,
        TraceRayParameters::InstanceMask,
        TraceRayParameters::HitGroup::Offset[RayType::Shadow],
        TraceRayParameters::HitGroup::GeometryStride,
        TraceRayParameters::MissShader::Offset[RayType::Shadow],
        rayDesc, shadowPayload);

    return shadowPayload.hit;
}

//***************************************************************************
//********************------ Ray gen shader.. -------************************
//***************************************************************************

uint getNewSeed(uint param1, uint param2, uint numPermutation)
{
    uint s0 = 0;
    uint v0 = param1;
    uint v1 = param2;

    for (uint perm = 0; perm < numPermutation; perm++)
    {
        s0 += 0x9e3779b9;
        v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s0) ^ ((v1 >> 5) + 0xc8013ea4);
        v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s0) ^ ((v0 >> 5) + 0x7e95761e);
    }

    return v0;
}

float rnd(float x)
{
    uint2 launchIdx = DispatchRaysIndex().xy;
    uint2 launchDim = DispatchRaysDimensions().xy;
    uint bufferOffset = launchDim.x * launchIdx.y + launchIdx.x;
    //uint seed = getNewSeed(bufferOffset, g_sceneCB.frameNum, 8);
    return frac(sin(bufferOffset * x + g_sceneCB.frameNum+ g_sceneCB.elapsedTime) * 43758.5453);
}
float3 random_in_unit_sphere()
{
    float3 a = { rnd(0.3),rnd(2.3),rnd(0.8) };
    float3 one = { 1,1,1 };
    return normalize(2 * a - one);
}
float3 sample_hemisphere_cos()
{
    float3 sampleDir;

    float param1 = rnd(0.8);
    float param2 = rnd(0.2);

    // Uniformly sample disk.
    float r = sqrt(param1);
    float phi = 2.0f * Pi * param2;
    sampleDir.x = r * cos(phi);
    sampleDir.y = r * sin(phi);

    // Project up to hemisphere.
    sampleDir.z = sqrt(max(0.0f, 1.0f - r * r));

    return sampleDir;
}

float3 sample_emitter()
{
    uint2 launchIdx = DispatchRaysIndex().xy;
    uint2 launchDim = DispatchRaysDimensions().xy;
    uint bufferOffset = launchDim.x * launchIdx.y + launchIdx.x;
    uint seed = getNewSeed(bufferOffset, bufferOffset, 8);

    return float3(1.5+2*rnd(1.3), 3 * rnd(1.2), 1.5 + 2 * rnd(11) );
}


[shader("raygeneration")]
void MyRaygenShader()
{
    // Generate a ray for a camera pixel corresponding to an index from the dispatched 2D grid.
 
    // Cast a ray into the scene and retrieve a shaded color.
    UINT currentRecursionDepth = 0;
    float4 color = { 0,0,0,0 };
	uint2 launchIdx = DispatchRaysIndex().xy;
    uint2 launchDim = DispatchRaysDimensions().xy;
    uint bufferOffset = launchDim.x * launchIdx.y + launchIdx.x;
    uint seed = getNewSeed(bufferOffset, bufferOffset, 8);
    [unroll]
    for (int i = 0; i < Sample_Num; i++) {
        float2 random = float2(rnd(0.9), rnd(5.4));
        Ray ray = GenerateCameraRay(DispatchRaysIndex().xy, g_sceneCB.cameraPosition.xyz, g_sceneCB.projectionToWorld);
        //GenerateCameraRay((float2)DispatchRaysIndex().xy , origin, rayDir);
        //Ray ray;
        //ray.origin = origin;
        //ray.direction = rayDir;
        color += TraceRadianceRay(ray, 0);
    }
    //float4 color = TraceRadianceRay(ray, currentRecursionDepth);

    // Write the raytraced color to the output texture.
    g_renderTarget[DispatchRaysIndex().xy] = color / Sample_Num;
}

//***************************************************************************
//******************------ Closest hit shaders -------***********************
//***************************************************************************

[shader("closesthit")]
void MyClosestHitShader_Triangle(inout RayPayload rayPayload, in BuiltInTriangleIntersectionAttributes attr)
{
    // Get the base index of the triangle's first 16 bit index.
    uint indexSizeInBytes = 2;
    uint indicesPerTriangle = 3;
    uint triangleIndexStride = indicesPerTriangle * indexSizeInBytes;
    uint baseIndex = PrimitiveIndex() * triangleIndexStride;

    // Load up three 16 bit indices for the triangle.
    const uint3 indices = Load3x16BitIndices(baseIndex, g_indices);

    // Retrieve corresponding vertex normals for the triangle vertices.
    float3 triangleNormal = g_vertices[indices[0]].normal;

    // PERFORMANCE TIP: it is recommended to avoid values carry over across TraceRay() calls. 
    // Therefore, in cases like retrieving HitWorldPosition(), it is recomputed every time.

    // Shadow component.
    // Trace a shadow ray.
    float3 hitPosition = HitWorldPosition();
    Ray shadowRay = { hitPosition, normalize(g_sceneCB.lightPosition.xyz - hitPosition) };
    bool shadowRayHit = TraceShadowRayAndReportIfHit(shadowRay, rayPayload.recursionDepth);

    float checkers = AnalyticalCheckersTexture(HitWorldPosition(), triangleNormal, g_sceneCB.cameraPosition.xyz, g_sceneCB.projectionToWorld);

    // Reflected component.
    float4 reflectedColor = float4(0, 0, 0, 0);
    float4 scatteredColor = float4(0, 0, 0, 0);
    float4 color = float4(0, 0, 0, 0);
    if (l_materialCB.reflectanceCoef > 0.001)
    {
        // Trace a reflection ray.
        Ray reflectionRay = { HitWorldPosition(), reflect(WorldRayDirection(), float3(0,1,0)) + l_materialCB.specularCoef * random_in_unit_sphere() };
        float4 reflectionColor = TraceRadianceRay(reflectionRay, rayPayload.recursionDepth);

        float3 fresnelR = FresnelReflectanceSchlick(WorldRayDirection(), float3(0, 1, 0), l_materialCB.albedo.xyz);
        reflectedColor = l_materialCB.reflectanceCoef * float4(fresnelR, 1) * reflectionColor;
        color = float4(reflectedColor.x * l_materialCB.albedo.x,
            reflectedColor.y * l_materialCB.albedo.y,
            reflectedColor.z * l_materialCB.albedo.z,
            reflectedColor.w * l_materialCB.albedo.w);
    }
    else if (l_materialCB.diffuseCoef > 0.001)
    {
        float light_rate = rnd(8.2);
        if (light_rate > 0.5) {
            float3 target = sample_emitter();
            Ray scatterRay = { HitWorldPosition(),target };
            scatteredColor = l_materialCB.diffuseCoef * TraceRadianceRay(scatterRay, rayPayload.recursionDepth);
            color = float4(scatteredColor.x * l_materialCB.albedo.x,
                scatteredColor.y * l_materialCB.albedo.y,
                scatteredColor.z * l_materialCB.albedo.z,
                scatteredColor.w * l_materialCB.albedo.w);
        }
        else {
            Ray scatterRay = { HitWorldPosition(),float3(0,1,0) + random_in_unit_sphere() };
            scatteredColor = l_materialCB.diffuseCoef * TraceRadianceRay(scatterRay, rayPayload.recursionDepth);
            color = float4(scatteredColor.x * l_materialCB.albedo.x,
                scatteredColor.y * l_materialCB.albedo.y,
                scatteredColor.z * l_materialCB.albedo.z,
                scatteredColor.w * l_materialCB.albedo.w);
        }
    }
    // Calculate final color.
    //float4 phongColor = CalculatePhongLighting(l_materialCB.albedo, attr.normal, shadowRayHit, l_materialCB.diffuseCoef, l_materialCB.specularCoef, l_materialCB.specularPower);

    //float4 color = l_materialCB.albedo + scatteredColor + reflectedColor;

    // Apply visibility falloff.
    float t = RayTCurrent();
    color = lerp(color, BackgroundColor, 1.0 - exp(-0.000002*t*t*t));

    rayPayload.color = color;
}

[shader("closesthit")]
void MyClosestHitShader_AABB(inout RayPayload rayPayload, in ProceduralPrimitiveAttributes attr)
{
    // PERFORMANCE TIP: it is recommended to minimize values carry over across TraceRay() calls. 
    // Therefore, in cases like retrieving HitWorldPosition(), it is recomputed every time.

    // Shadow component.
    // Trace a shadow ray.
    float3 hitPosition = HitWorldPosition();
    Ray shadowRay = { hitPosition, normalize(g_sceneCB.lightPosition.xyz - hitPosition) };
    bool shadowRayHit = TraceShadowRayAndReportIfHit(shadowRay, rayPayload.recursionDepth);

    // Reflected component.
    float4 reflectedColor = float4(0, 0, 0, 0);
    float4 scatteredColor = float4(0, 0, 0, 0);
    float4 color = float4(0, 0, 0, 0);
    if(l_materialCB.emission.x>0.1){
        rayPayload.color = l_materialCB.emission;
        return;
    }
    if (l_materialCB.reflectanceCoef > 0.001)
    {
        // Trace a reflection ray.
        Ray reflectionRay = { HitWorldPosition(), reflect(WorldRayDirection(), attr.normal) + l_materialCB.specularCoef * random_in_unit_sphere()};
        float4 reflectionColor = TraceRadianceRay(reflectionRay, rayPayload.recursionDepth);

        float3 fresnelR = FresnelReflectanceSchlick(WorldRayDirection(), attr.normal, l_materialCB.albedo.xyz);
        reflectedColor = l_materialCB.reflectanceCoef * float4(fresnelR, 1) * reflectionColor;
        color = reflectedColor;
    }
    else if(l_materialCB.diffuseCoef > 0.001)
    {
        float light_rate = rnd(123.2);
        if (light_rate > 0.5) {
            float3 target = sample_emitter();
            Ray scatterRay = { HitWorldPosition(), target };
            scatteredColor = l_materialCB.diffuseCoef * TraceRadianceRay(scatterRay, rayPayload.recursionDepth);
            color = float4(scatteredColor.x * l_materialCB.albedo.x,
                scatteredColor.y * l_materialCB.albedo.y,
                scatteredColor.z * l_materialCB.albedo.z,
                scatteredColor.w * l_materialCB.albedo.w);
        }
        else {
            Ray scatterRay = { HitWorldPosition(),attr.normal + random_in_unit_sphere() };
            scatteredColor = l_materialCB.diffuseCoef * TraceRadianceRay(scatterRay, rayPayload.recursionDepth);
            color = float4(scatteredColor.x * l_materialCB.albedo.x,
                scatteredColor.y * l_materialCB.albedo.y,
                scatteredColor.z * l_materialCB.albedo.z,
                scatteredColor.w * l_materialCB.albedo.w);
        }
    }
    // Calculate final color.
    //float4 phongColor = CalculatePhongLighting(l_materialCB.albedo, attr.normal, shadowRayHit, l_materialCB.diffuseCoef, l_materialCB.specularCoef, l_materialCB.specularPower);
    
    //color = l_materialCB.albedo + scatteredColor + reflectedColor;

    // Apply visibility falloff.
    float t = RayTCurrent();
    color = lerp(color, BackgroundColor, 1.0 - exp(-0.000002*t*t*t));

    rayPayload.color = color + l_materialCB.emission;
}

//***************************************************************************
//**********************------ Miss shaders -------**************************
//***************************************************************************

[shader("miss")]
void MyMissShader(inout RayPayload rayPayload)
{
    float4 backgroundColor = float4(BackgroundColor);
    rayPayload.color = backgroundColor;
}

[shader("miss")]
void MyMissShader_ShadowRay(inout ShadowRayPayload rayPayload)
{
    rayPayload.hit = false;
}

//***************************************************************************
//*****************------ Intersection shaders-------************************
//***************************************************************************

// Get ray in AABB's local space.
Ray GetRayInAABBPrimitiveLocalSpace()
{
    PrimitiveInstancePerFrameBuffer attr = g_AABBPrimitiveAttributes[l_aabbCB.instanceIndex];

    // Retrieve a ray origin position and direction in bottom level AS space 
    // and transform them into the AABB primitive's local space.
    Ray ray;
    ray.origin = mul(float4(ObjectRayOrigin(), 1), attr.bottomLevelASToLocalSpace).xyz;
    ray.direction = mul(ObjectRayDirection(), (float3x3) attr.bottomLevelASToLocalSpace);
    return ray;
}

[shader("intersection")]
void MyIntersectionShader_AnalyticPrimitive()
{
    Ray localRay = GetRayInAABBPrimitiveLocalSpace();
    AnalyticPrimitive::Enum primitiveType = (AnalyticPrimitive::Enum) l_aabbCB.primitiveType;

    float thit;
    ProceduralPrimitiveAttributes attr = (ProceduralPrimitiveAttributes)0;
    if (RayAnalyticGeometryIntersectionTest(localRay, primitiveType, thit, attr))
    {
        PrimitiveInstancePerFrameBuffer aabbAttribute = g_AABBPrimitiveAttributes[l_aabbCB.instanceIndex];
        attr.normal = mul(attr.normal, (float3x3) aabbAttribute.localSpaceToBottomLevelAS);
        attr.normal = normalize(mul((float3x3) ObjectToWorld3x4(), attr.normal));

        ReportHit(thit, /*hitKind*/ 0, attr);
    }
}

[shader("intersection")]
void MyIntersectionShader_VolumetricPrimitive()
{
    Ray localRay = GetRayInAABBPrimitiveLocalSpace();
    VolumetricPrimitive::Enum primitiveType = (VolumetricPrimitive::Enum) l_aabbCB.primitiveType;
    
    float thit;
    ProceduralPrimitiveAttributes attr = (ProceduralPrimitiveAttributes)0;
    if (RayVolumetricGeometryIntersectionTest(localRay, primitiveType, thit, attr, g_sceneCB.elapsedTime))
    {
        PrimitiveInstancePerFrameBuffer aabbAttribute = g_AABBPrimitiveAttributes[l_aabbCB.instanceIndex];
        attr.normal = mul(attr.normal, (float3x3) aabbAttribute.localSpaceToBottomLevelAS);
        attr.normal = normalize(mul((float3x3) ObjectToWorld3x4(), attr.normal));

        ReportHit(thit, /*hitKind*/ 0, attr);
    }
}

[shader("intersection")]
void MyIntersectionShader_SignedDistancePrimitive()
{
    Ray localRay = GetRayInAABBPrimitiveLocalSpace();
    SignedDistancePrimitive::Enum primitiveType = (SignedDistancePrimitive::Enum) l_aabbCB.primitiveType;

    float thit;
    ProceduralPrimitiveAttributes attr = (ProceduralPrimitiveAttributes)0;
    if (RaySignedDistancePrimitiveTest(localRay, primitiveType, thit, attr, l_materialCB.stepScale))
    {
        PrimitiveInstancePerFrameBuffer aabbAttribute = g_AABBPrimitiveAttributes[l_aabbCB.instanceIndex];
        attr.normal = mul(attr.normal, (float3x3) aabbAttribute.localSpaceToBottomLevelAS);
        attr.normal = normalize(mul((float3x3) ObjectToWorld3x4(), attr.normal));
        
        ReportHit(thit, /*hitKind*/ 0, attr);
    }
}

#endif // RAYTRACING_HLSL