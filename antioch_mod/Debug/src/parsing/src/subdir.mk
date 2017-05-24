################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../src/parsing/src/ascii_parser.C \
../src/parsing/src/blottner_parsing.C \
../src/parsing/src/cea_mixture_ascii_parsing.C \
../src/parsing/src/cea_mixture_parsing.C \
../src/parsing/src/chemkin_parser.C \
../src/parsing/src/nasa_mixture_parsing.C \
../src/parsing/src/parser_base.C \
../src/parsing/src/read_reaction_set_data.C \
../src/parsing/src/species_parsing.C \
../src/parsing/src/sutherland_parsing.C \
../src/parsing/src/transport_species_parsing.C \
../src/parsing/src/xml_parser.C 

O_SRCS += \
../src/parsing/src/ascii_parser.o \
../src/parsing/src/blottner_parsing.o \
../src/parsing/src/cea_mixture_ascii_parsing.o \
../src/parsing/src/cea_mixture_parsing.o \
../src/parsing/src/chemkin_parser.o \
../src/parsing/src/nasa_mixture_parsing.o \
../src/parsing/src/parser_base.o \
../src/parsing/src/read_reaction_set_data.o \
../src/parsing/src/species_parsing.o \
../src/parsing/src/sutherland_parsing.o \
../src/parsing/src/transport_species_parsing.o \
../src/parsing/src/xml_parser.o 

C_UPPER_DEPS += \
./src/parsing/src/ascii_parser.d \
./src/parsing/src/blottner_parsing.d \
./src/parsing/src/cea_mixture_ascii_parsing.d \
./src/parsing/src/cea_mixture_parsing.d \
./src/parsing/src/chemkin_parser.d \
./src/parsing/src/nasa_mixture_parsing.d \
./src/parsing/src/parser_base.d \
./src/parsing/src/read_reaction_set_data.d \
./src/parsing/src/species_parsing.d \
./src/parsing/src/sutherland_parsing.d \
./src/parsing/src/transport_species_parsing.d \
./src/parsing/src/xml_parser.d 

OBJS += \
./src/parsing/src/ascii_parser.o \
./src/parsing/src/blottner_parsing.o \
./src/parsing/src/cea_mixture_ascii_parsing.o \
./src/parsing/src/cea_mixture_parsing.o \
./src/parsing/src/chemkin_parser.o \
./src/parsing/src/nasa_mixture_parsing.o \
./src/parsing/src/parser_base.o \
./src/parsing/src/read_reaction_set_data.o \
./src/parsing/src/species_parsing.o \
./src/parsing/src/sutherland_parsing.o \
./src/parsing/src/transport_species_parsing.o \
./src/parsing/src/xml_parser.o 


# Each subdirectory must supply rules for building sources it contributes
src/parsing/src/%.o: ../src/parsing/src/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../ -I/home/bbojko/local/mpich-3.2/include -I/home/bbojko/Projects/GrossDiffFlame/petsc-2.2.1/include -I/home/bbojko/Projects/GrossDiffFlame/petsc-2.2.1/bmake/linux-gnu -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/utilities/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/units/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/transport/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/kinetics/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/thermo/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/thermal_conduction/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/particles_flux/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/diffusion/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/parsing/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/core/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/viscosity/include" -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


