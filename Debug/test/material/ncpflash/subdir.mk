################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/material/ncpflash/test_ncpflash.cc 

CC_DEPS += \
./test/material/ncpflash/test_ncpflash.d 

OBJS += \
./test/material/ncpflash/test_ncpflash.o 


# Each subdirectory must supply rules for building sources it contributes
test/material/ncpflash/%.o: ../test/material/ncpflash/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


