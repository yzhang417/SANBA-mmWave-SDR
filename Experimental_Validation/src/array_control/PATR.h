#ifdef MATLAB_MEX_FILE
#include <tmwtypes.h>
#else
#include "rtwtypes.h"
#endif

int OpenSerialPortRx(const char *portName);

int CloseSerialPortRx(void);

int WriteCommandRx(const char *input_command, const char *portName);

int ReadCommandRx(void); 

int ReadCommandRxPrint(void);

int OpenSerialPortTx(const char *portName);

int CloseSerialPortTx(void);

int WriteCommandTx(const char *input_command, const char *portName);

int ReadCommandTx(void);

int ReadCommandTxPrint(void);

int Operation_Interval_Time(int nMicroSecond);