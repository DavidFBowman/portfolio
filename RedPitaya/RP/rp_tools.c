/**
	Set of functions to drive the Red Pitaya signal acquisition board. 
	The base of these functions are provided by Red Pitaya and I have rewritten them to compile as a python library.

*/
#include <python2.7/Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "./rp.h"
#include "./rp_tools.h"


/**
This MUST be run before any other functions. 
*/
void  Initialize ()
{ 
    rp_Init();
}


/**
Run this when you have finished acquisition.
*/
void Finalize ()
{
    rp_Reset();
}

/**
Read out from the buffer in port IN1.
Change the value in the sleep function to be compatible with the time taken to fill the buffer.
1 second is fine as a first estimate. 
*/
float *acquire_buffer(){

        uint32_t buff_size = 16384;
        float *buff = (float *)malloc(buff_size * sizeof(float));

        rp_AcqReset();
        rp_AcqSetDecimation(1);
        rp_AcqSetTriggerDelay(0);

        rp_AcqStart();


        sleep(1);
        rp_AcqSetTriggerSrc(RP_TRIG_SRC_NOW);
        rp_acq_trig_state_t state = RP_TRIG_STATE_TRIGGERED;

        while(1){
                rp_AcqGetTriggerState(&state);
                if(state == RP_TRIG_STATE_TRIGGERED){
                sleep(1);
                break;
                }
        }

        rp_AcqGetOldestDataV(RP_CH_1, &buff_size, buff);


        return buff;

}

/**
Generate a sine wave locally from port OUT1.
*/
void sine_wave(){

    /* Print error, if rp_Init() function failed */
    static int i;
    if (i == 0)
    {   
        ++i;
        Initialize();
        atexit(Finalize);
    }



    /* Generating frequency */
    rp_GenFreq(RP_CH_1, 10000.0);

    /* Generating amplitude */
    rp_GenAmp(RP_CH_1, 1.0);

    /* Generating wave form */
    rp_GenWaveform(RP_CH_1, RP_WAVEFORM_SINE);

    /* Enable channel */
    rp_GenOutEnable(RP_CH_1);


}

/**
Kill the output from port OUT1
*/
void stop_signal(){


	rp_GenOutDisable(RP_CH_1);


}