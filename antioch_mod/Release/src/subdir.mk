################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/antioch_mod.cpp 

OBJS += \
./src/antioch_mod.o 

CPP_DEPS += \
./src/antioch_mod.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../ -I/home/bbojko/local/mpich-3.2/include -I/home/bbojko/Projects/GrossDiffFlame/petsc-2.2.1/include -I/home/bbojko/Projects/GrossDiffFlame/petsc-2.2.1/bmake/linux-gnu -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/core/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/diffusion/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/kinetics/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/parsing/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/particles_flux/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/thermal_conduction/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/thermo/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/transport/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/units/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/utilities/include" -I"/home/bbojko/Workspaces/antioch_mod/antioch_mod/src/viscosity/include" -O3 -pg -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


