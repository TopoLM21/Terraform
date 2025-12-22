if (NOT DEFINED QWT_SOURCE_DIR)
    message(FATAL_ERROR "QWT_SOURCE_DIR is required.")
endif()

if (NOT DEFINED QWT_INSTALL_PREFIX)
    message(FATAL_ERROR "QWT_INSTALL_PREFIX is required.")
endif()

set(qwtconfig_path "${QWT_SOURCE_DIR}/qwtconfig.pri")

if (NOT EXISTS "${qwtconfig_path}")
    message(FATAL_ERROR "Qwt config file not found at ${qwtconfig_path}.")
endif()

file(READ "${qwtconfig_path}" qwtconfig_contents)

string(REGEX REPLACE
    "^[ \t]*QWT_INSTALL_PREFIX[ \t]*=.*$"
    "QWT_INSTALL_PREFIX = ${QWT_INSTALL_PREFIX}"
    qwtconfig_updated
    "${qwtconfig_contents}"
)

if (qwtconfig_updated STREQUAL qwtconfig_contents)
    message(FATAL_ERROR "QWT_INSTALL_PREFIX entry was not found in ${qwtconfig_path}.")
endif()

file(WRITE "${qwtconfig_path}" "${qwtconfig_updated}")
