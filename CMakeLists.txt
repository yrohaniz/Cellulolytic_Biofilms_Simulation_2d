cmake_minimum_required(VERSION 3.8)
project(Cellulolytic_Biofilms_Simulation)

set(CMAKE_C_STANDARD 11)

set(SOURCE_FILES
        main.c
        model_functions.c
        CeBiSi.h
        integrator.c
        approximations.c
        PRNG_uniform.c
        Marsaglia_polar_generator.c
        LCGS.c
        params_read.c pararand.c)
add_executable(Cellulolytic_Biofilms_Simulation ${SOURCE_FILES})

target_link_libraries(Cellulolytic_Biofilms_Simulation m)