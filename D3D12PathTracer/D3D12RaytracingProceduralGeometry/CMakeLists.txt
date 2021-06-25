set(PROJECT_NAME D3D12RaytracingProceduralGeometry)

################################################################################
# Source groups
################################################################################
set(no_group_source_files
    "readme.md"
)
source_group("" FILES ${no_group_source_files})

set(Assets__Shaders
    "RayTracingHlslCompat.h"
    "RaytracingShaderHelper.hlsli"
)
source_group("Assets\\Shaders" FILES ${Assets__Shaders})

set(Assets__Shaders__ProceduralPrimitives
    "AnalyticPrimitives.hlsli"
    "ProceduralPrimitivesLibrary.hlsli"
    "SignedDistanceFractals.hlsli"
    "SignedDistancePrimitives.hlsli"
    "VolumetricPrimitives.hlsli"
)
source_group("Assets\\Shaders\\ProceduralPrimitives" FILES ${Assets__Shaders__ProceduralPrimitives})

set(Assets__Shaders__Util
    "util/HlslCompat.h"
)
source_group("Assets\\Shaders\\Util" FILES ${Assets__Shaders__Util})

set(Header_Files
    "D3D12RaytracingProceduralGeometry.h"
    "DirectXRaytracingHelper.h"
    "RaytracingSceneDefines.h"
    "stdafx.h"
)
source_group("Header Files" FILES ${Header_Files})

set(Header_Files__Util
    "util/DeviceResources.h"
    "util/DXSample.h"
    "util/DXSampleHelper.h"
    "util/PerformanceTimers.h"
    "util/StepTimer.h"
    "util/Win32Application.h"
)
source_group("Header Files\\Util" FILES ${Header_Files__Util})

set(Source_Files
    "D3D12RaytracingProceduralGeometry.cpp"
    "Main.cpp"
    "stdafx.cpp"
)
source_group("Source Files" FILES ${Source_Files})

set(Source_Files__Util
    "util/DeviceResources.cpp"
    "util/DXSample.cpp"
    "util/PerformanceTimers.cpp"
    "util/Win32Application.cpp"
)
source_group("Source Files\\Util" FILES ${Source_Files__Util})

set(ALL_FILES
    ${no_group_source_files}
    ${Assets__Shaders}
    ${Assets__Shaders__ProceduralPrimitives}
    ${Assets__Shaders__Util}
    ${Header_Files}
    ${Header_Files__Util}
    ${Source_Files}
    ${Source_Files__Util}
)

################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${ALL_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Samples/Tutorials")

add_precompiled_header(${PROJECT_NAME} "stdafx.h" "stdafx.cpp")

use_props(${PROJECT_NAME} "${CMAKE_CONFIGURATION_TYPES}" "${DEFAULT_CXX_PROPS}")
set(ROOT_NAMESPACE D3D12Raytracing)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VS_GLOBAL_KEYWORD "Win32Proj"
)
################################################################################
# Output directory
################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES
    OUTPUT_DIRECTORY_DEBUG   "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/"
    OUTPUT_DIRECTORY_PROFILE "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/"
    OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_VS_PLATFORM_NAME}/$<CONFIG>/"
)
set_target_properties(${PROJECT_NAME} PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION_PROFILE "TRUE"
    INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
)
################################################################################
# Include directories
################################################################################
target_include_directories(${PROJECT_NAME} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/util;"
    "${CMAKE_CURRENT_SOURCE_DIR}/..;"
    "${CMAKE_CURRENT_BINARY_DIR}\${CMAKE_CFG_INTDIR}\"
)

################################################################################
# Compile definitions
################################################################################
target_compile_definitions(${PROJECT_NAME} PRIVATE
    "$<$<CONFIG:Debug>:"
        "_DEBUG"
    ">"
    "$<$<CONFIG:Profile>:"
        "NDEBUG"
    ">"
    "$<$<CONFIG:Release>:"
        "NDEBUG"
    ">"
    "WIN32;"
    "_WINDOWS;"
    "UNICODE;"
    "_UNICODE"
)

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE
        $<$<CONFIG:Debug>:
            /Od
        >
        $<$<CONFIG:Profile>:
            /O2;
            /Oi;
            /Gy
        >
        $<$<CONFIG:Release>:
            /O2;
            /Oi;
            /Gy
        >
        /sdl;
        /W3;
        ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
        ${DEFAULT_CXX_EXCEPTION_HANDLING}
    )
    target_link_options(${PROJECT_NAME} PRIVATE
        $<$<CONFIG:Debug>:
            /INCREMENTAL
        >
        $<$<CONFIG:Profile>:
            /OPT:REF;
            /OPT:ICF;
            /INCREMENTAL:NO
        >
        $<$<CONFIG:Release>:
            /OPT:REF;
            /OPT:ICF;
            /INCREMENTAL:NO
        >
        /DELAYLOAD:d3d12.dll;
        /DEBUG;
        /SUBSYSTEM:WINDOWS
    )
endif()

################################################################################
# Dependencies
################################################################################
set(ADDITIONAL_LIBRARY_DEPENDENCIES
    "$<$<CONFIG:Debug>:"
        "dxguid"
    ">"
    "d3d12;"
    "dxgi;"
    "kernel32;"
    "user32;"
    "gdi32;"
    "winspool;"
    "comdlg32;"
    "advapi32;"
    "shell32;"
    "ole32;"
    "oleaut32;"
    "uuid;"
    "odbc32;"
    "odbccp32"
)
target_link_libraries(${PROJECT_NAME} PRIVATE "${ADDITIONAL_LIBRARY_DEPENDENCIES}")

use_package(${PROJECT_NAME} WinPixEventRuntime 1.0.180612001)
