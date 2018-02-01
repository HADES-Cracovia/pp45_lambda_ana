# - Find HYDRA instalation
# F.Uhlig@gsi.de (fairroot.gsi.de)

# Oryginal module written for ROOT, adapted for HYDRA by R.Lalik@ph.tum.de


MESSAGE(STATUS "Looking for Hydra...")

IF (NOT HADDIR)
    SET (INTHADDIR $ENV{HADDIR})
ENDIF (NOT HADDIR)

IF (NOT INTHADDIR)
MESSAGE( FATAL_ERROR "HADDIR is not set. Please set HADDIR first.")
ENDIF (NOT INTHADDIR)

IF (NOT MYHADDIR)
        SET (MYHADDIR $ENV{MYHADDIR})
ENDIF (NOT MYHADDIR)

SET(HYDRA_PACKAGE_SEARCHPATH
    ${INTHADDIR}/lib
)

SET(HYDRA_DEFINITIONS "")

SET(HYDRA_INSTALLED_VERSION_TOO_OLD FALSE)

SET(HYDRA_MAIN_LIBRARY HYDRA_MAIN_LIBRARY-NOTFOUND)

FIND_LIBRARY(HYDRA_MAIN_LIBRARY NAMES Pid PATHS
    ${HYDRA_PACKAGE_SEARCHPATH}
    NO_DEFAULT_PATH)
    
IF (${HYDRA_MAIN_LIBRARY} MATCHES "HYDRA_MAIN_LIBRARY-NOTFOUND")
    MESSAGE( STATUS "Hydra library not found. Please check your Hydra installation.")
    SET(HYDRA_FOUND FALSE)
ELSE (${HYDRA_MAIN_LIBRARY} MATCHES "HYDRA_MAIN_LIBRARY-NOTFOUND")
    MESSAGE(STATUS "Looking for Hydra... - found ${INTHADDIR}")
    MESSAGE(STATUS "   MYHADDIR = ${MYHADDIR}")
    SET(HYDRA_FOUND TRUE)
ENDIF (${HYDRA_MAIN_LIBRARY} MATCHES "HYDRA_MAIN_LIBRARY-NOTFOUND")  

IF (HYDRA_FOUND)
    SET(HYDRA_LIBRARY_DIR ${INTHADDIR}/lib )
    SET(HYDRA_INCLUDE_DIR ${INTHADDIR}/include )

    set(HYDRA_LIBRARIES)
    foreach(_cpt Alignment Dst HadesGo4 Hodo Hydra Hyp Kick MdcGarfield MdcPid Mdc MdcTrackD MdcTrackG MdcTrackS MdcUtil Pairs PhyAna Pid PidUtil QA Revt Rich RichUtil Rpc Shower ShowerTofino ShowerUtil Simulation Start Tofino Tof TofUtil Tools Trigger TriggerUtil Wall Ora OraSim OraUtil)
        find_library(HYDRA_${_cpt}_LIBRARY NAMES ${_cpt} lib{_cpt} HINTS ${HYDRA_LIBRARY_DIR})
        if(HYDRA_${_cpt}_LIBRARY)
            mark_as_advanced(HYDRA_${_cpt}_LIBRARY)
            list(APPEND HYDRA_LIBRARIES ${HYDRA_${_cpt}_LIBRARY})
        endif()
    endforeach()
    if(HYDRA_LIBRARIES)
        list(REMOVE_DUPLICATES HYDRA_LIBRARIES)
    endif()

    # Make variables changeble to the advanced user
    MARK_AS_ADVANCED( HYDRA_LIBRARY_DIR HYDRA_INCLUDE_DIR HYDRA_DEFINITIONS)

    # Set HYDRA_INCLUDES
    SET(HYDRA_INCLUDES ${HYDRA_INCLUDE_DIR})

    SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${HYDRA_LIBRARY_DIR})
ENDIF (HYDRA_FOUND)
