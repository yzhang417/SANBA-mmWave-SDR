## Explanation of data set

- cal_data/
  This folder contains the collected data which will be used to
  calculate the radiation gain and phase error of each antenna 
  element. 

      Data in cal/data/Element_gain is collected by 
      [Experimental_Validation/main_programs/general_receiver_decoder/
      /sdruQPSKReceiver.m] with Program_ID = 9;

      Data in cal/data/Noise_floor is collected by 
      [Experimental_Validation/main_programs/fine_phased_array_calibration
      /estimate_noise_floor.m] and removing [clear all] command in 
      [Experimental_Validation/main_programs/general_receiver_decoder/
      /sdruQPSKReceiver.m].

      Data in cal/data/Tx_cal is collected by 
      [Experimental_Validation/main_programs/general_receiver_decoder/
      /sdruQPSKReceiver.m] with Program_ID = 9.1;

      Data in cal/data/Rx_Cal is collected by 
      [Experimental_Validation/main_programs/general_receiver_decoder/
      /sdruQPSKReceiver.m] with Program_ID = 9.2;

- cal_result/ 
  This folder contains the calibration result after running the 
  functions in [Experimental_Validation/main_programs/
  fine_phased_array_calibration/*.m].

- val_cal_data/ 
  This folder contains the collected data which will be used to
  validate the proposed calibration method. The data is collected
  by [Experimental_Validation/main_programs/main_cal_val_receiver.m].

- val_cpr_data/ 
  This folder contains the collected data which will be used to
  validate the proposed two-stage non-coherent beam alignment.
  The data is collected by [Experimental_Validation/main_programs/
  main_two_stage_cpr_receiver.m].

