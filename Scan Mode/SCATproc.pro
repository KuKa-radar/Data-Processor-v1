@SCATlib
PRO SCATproc                    ;this is the main program
common waitkeyvars
;
;rewriting SCATproc for Ku/Ka-SCATs
;new features: 
; processes either scanning or non-scanning files 
;author: J. Mead
;date: 6/13/19
;updated July 2020 to add batch processing option and processing limits on
;azimuth and elevation axis
;
  device,decompose=0
  NEXRAD_ct

  while !D.window ne -1 do wdelete
;select file and read file header
  !EXCEPT=2                     ;report all math errors as they occur
read,'enter 1 for Ku-Scat; 2 for Ka-Scat: ',iscat

; control variables
  read_configuration,configvars,iscat
  instrument_name=configvars.instrument_name
  calfile=configvars.calfile    ;pass this as plain variable since it may get modified if creating a new calfile
  raw_data_path=configvars.raw_data_path
  processed_data_path=configvars.processed_data_path
  if configvars.i_batch eq 0 then begin
     filename=dialog_pickfile(path=raw_data_path,/read,filter="*dat",title="choose SCAT raw data file")
     filecount=1
     files=strarr(1)
     files(0)=filename
  endif else begin
     files = file_search(raw_data_path,'*.dat',count=filecount) 
  endelse

  for ifile=0,filecount-1 do begin
     filename = files(ifile)
     print,filename
     readheader,filename,raw_data_path,scatvars,calvars,private_config,file_header_size
     print,'scatvars:',scatvars
     print,'calvars:',calvars
     print,'private_config:',private_config
     print,'file_header_size:',file_header_size
     
;
;read raw data file
;

     readraw,configvars,raw,scan_index,sweep_count,transition_flag,elevation,n_blocks,n_pol,scatvars,filename,file_header_size,$
             elapsed_time,time_sec,gps_latitude,gps_longitude,along_track_tilt,cross_track_tilt,independent_sample_index,distance,az_proc_index,sweep_count_override
;scan_index(*)=-4
     startchar=19
     filename_length = 16
     base_filename=instrument_name+strmid(filename,startchar,filename_length,/reverse_offset)
     processed_data_filename = processed_data_path+'/'+base_filename+'.nrcs'
     print,'processed_data_filename: ',processed_data_filename
     

;
;compute range profiles of reflectivity from raw data
;

     compute_range_profiles,range_gate_spacing,n_blocks,n_pol,ngates,gate_offset,gate_plot_max,elevation,scan_index,$
                            sweep_count,raw,pos_height,height,range,spec,gate_max_cal,geometric_peak,scatvars,private_config,configvars

;the calibrate procedure is used to generate the transmit distortion
;matrix 

     calibrate,calvars,configvars,gate_offset,gate_plot_max,spec,ngates,n_blocks,n_pol,range,gate_max_cal,reference_calibration_loop_power,$
               current_calibration_loop_power,ainv,finv,range_gate_spacing,corner_range,scatvars,scan_index,sweep_count,$
               elevation,total_corner_power_vv,total_corner_power_hh,base_filename,calfile

;
;compute calibrated covariance matrix and Mueller matrix containing average polarimetric
;scattering properties of target
;

     if max(scan_index) eq 1 then  polarimetric_processing_scan,scatvars,configvars,calvars,geometric_peak,range_peak_signal,scan_index,$
        sweep_count,transition_flag,elevation,line_elevation,line_height,ngates,nlines,$
        range_gate_spacing,range,height,gate_offset,gate_plot_max,spec,ainv,finv,l_matrix,$
        c_matrix,l_chi_gamma,positioner_state,rho_hv_vec,phase_hv_deg_vec,total_power,range_centroid_signal,az_proc_index,sweep_count_override

     if max(scan_index) ne 1 and configvars.i_proc_ind eq 0 then   polarimetric_processing_stare,scatvars,configvars,calvars,ngates,range_gate_spacing,$
        range,gate_offset,gate_plot_max,spec,ainv,finv,l_matrix,c_matrix,rho_hv_vec,phase_hv_deg_vec,total_power,n_blocks,$
        elapsed_time,n_groups,range_peak_signal,range_centroid_signal,$
        group_index,n_blocks_per_group,scan_index
     if max(scan_index) ne 1 and configvars.i_proc_ind eq 1 then polarimetric_processing_stare_by_independent_samples,scatvars,$
        configvars,calvars,ngates,range_gate_spacing,range,gate_offset,gate_plot_max,spec,ainv,finv,l_matrix,c_matrix,$
        rho_hv_vec,phase_hv_deg_vec,total_power,n_blocks,elapsed_time,n_groups,range_peak_signal,range_centroid_signal,$
        group_index,n_blocks_per_group,scan_index,independent_sample_index,distance 

                                ;
;compute nrcs and store processed data. Computes NRCS from covariance
;matrix elements
;
;nrcs_compute uses formula from T. Geldsetzer, J.B. Mead,
;et al., "Surface-Based Polarimetric C-band Scatterometer for
;Field Measurements of Sea Ice", IEEE Transactions on
;Geoscience and Remote Sensing, November 2007, pg. 3405-3415, equation (1)  
;

     if max(scan_index) eq 1 then nrcs_compute_scan,scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,$
        range_centroid_signal,corner_range,pos_height,line_elevation,nlines,reference_calibration_loop_power,$
        current_calibration_loop_power,L_matrix,processed_data_filename,rho_hv_vec,$
        phase_hv_deg_vec,total_power,line_height,total_corner_power_vv,total_corner_power_hh,gps_latitude,gps_longitude,calfile

     if max(scan_index) ne 1 then nrcs_compute_stare,scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,$
        range_centroid_signal,corner_range,pos_height,n_groups,reference_calibration_loop_power,current_calibration_loop_power,$
        L_matrix,processed_data_filename,rho_hv_vec,phase_hv_deg_vec,total_power,total_corner_power_vv,total_corner_power_hh,$
        calfile,time_sec,gps_latitude,gps_longitude,group_index,along_track_tilt,cross_track_tilt,n_blocks_per_group,independent_sample_index

  endfor
  
END 
SCATproc;executes procedure scatproc
end


