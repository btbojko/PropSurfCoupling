################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../src/utilities/src/antioch_version.C \
../src/utilities/src/gsl_spliner_impl.C \
../src/utilities/src/gsl_spliner_shim.C \
../src/utilities/src/string_utils.C 

O_SRCS += \
../src/utilities/src/antioch_version.o \
../src/utilities/src/gsl_spliner_impl.o \
../src/utilities/src/gsl_spliner_shim.o \
../src/utilities/src/string_utils.o 

C_UPPER_DEPS += \
./src/utilities/src/antioch_version.d \
./src/utilities/src/gsl_spliner_impl.d \
./src/utilities/src/gsl_spliner_shim.d \
./src/utilities/src/string_utils.d 

OBJS += \
./src/utilities/src/antioch_version.o \
./src/utilities/src/gsl_spliner_impl.o \
./src/utilities/src/gsl_spliner_shim.o \
./src/utilities/src/string_utils.o 


# Each subdirectory must supply rules for building sources it contributes
src/utilities/src/%.o: ../src/utilities/src/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../ -I/home/bbojko/local/mpich-3.2/include -I/home/bbojko/Projects/GrossDiffFlame/petsc-2.2.1/include -I/home/bbojko/Projects/GrossDiffFlame/petsc-2.2.1/bmake/linux-gnu -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/utilities/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/units/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/transport/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/kinetics/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/thermo/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/thermal_conduction/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/particles_flux/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/diffusion/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/parsing/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/core/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/viscosity/include" -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


