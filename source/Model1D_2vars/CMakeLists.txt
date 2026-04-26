# Definujeme jméno binárky pro tento model (např. 'run_model_1d')
set(TARGET "model_1d")

# Seznam souborů (cesty jsou relativní k TOMUTO CMakeLists)
set(TARGET_SRC
    model.cc
    ../../main.cc  # main.cc je o dvě úrovně výš
)

# Vytvoření binárky
add_executable(${TARGET} ${TARGET_SRC})

# Nastavení deal.II pro tuto binárku
deal_ii_setup_target(${TARGET})

# Aby kód viděl na AceGen složku uvnitř tohoto modelu
target_include_directories(${TARGET} PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/AceGen
)

# Kopírování parametrů do build/source/Model1D_2vars/
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/double_ditch.prm
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})