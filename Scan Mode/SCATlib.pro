Pro waitkey,dum
common waitkeyvars,configvars
;uncomment the following two lines when compiling SCATproc.sav for
;use when running unlicensed IDL
; SPW_CONTINUE
if configvars.i_batch eq 1 then goto,jump_waitkey
  print,'press any key to continue, q to quit'
  dum=string(0)
  dum = get_kbrd(1)
  IF dum EQ 'q' THEN stop
  IF dUM EQ 'Q' THEN stop
jump_waitkey:
END;waitkey

PRO SPW_CONTINUE
;***create base***
cont_base = WIDGET_BASE(TITLE='Click to Continue:', /FRAME)
;***create button:***
continue_but= WIDGET_BUTTON(cont_base, /FRAME,$
VALUE=' CONTINUE ',UVALUE='CONT_CONT')
;***Realize the menu
WIDGET_CONTROL, cont_base, /REALIZE
;*** Wait for the first event: ***
event=WIDGET_EVENT(cont_base)

;*** Destroy the GUI, let caller program proceed: ***
WIDGET_CONTROL, cont_base, /DESTROY

END; SPW_CONTINUE

function counts_to_voltage,scatvars,counts,i_analyze_smoothing 
voltage=(counts/51.4)/1000.;volts
return,voltage
end

PRO readheader,filename,raw_data_path,scatvars,calvars,private_config,file_header_size

@file_header.txt


fd = 5
close,fd
openu,fd,filename
readu,fd,file_header_t


scatvars=file_header_t.config
calvars=file_header_t.calibration
private_config=file_header_t.private_config
print
print,'AD clock frequency (MHz): ',scatvars.ad_clock_frequency_hz/1e06
print,'AD trig delay (ns): ',private_config.ad_trig_delay_ns
print,'chirp bandwidth (MHz): ',scatvars.chirp_bandwidth_hz/1e06 
print,'chirp center frequency (MHz): ',scatvars.chirp_center_frequency_hz/1e06
print,'chirp amplitude dBm: ',scatvars.chirp_amplitude_dbm
print,'chirp width (ms): ',scatvars.chirp_width_ns/1e06
print,'decimation: ',scatvars.decimation
print,'decimation mode: ',scatvars.decimation_mode
print,'frame delay (ms): ',scatvars.frame_delay_ns/1e06

;
;in the context of the server code, n_gates is the number of samples
;per sweep. In the IDL code, after the FFT, I use the variable ngates
;(without underbar) to denote the number of range gates
;

print,'number of time samples per sweep (n_gates): ',scatvars.n_gates
print,'pulse repetition period (ms): ',scatvars.prp_ns/1e06
print,'nominal range resolution (c_light/2/chirp_bandwidth, m): ',scatvars.range_resolution_m
print,'server state, 0=idle; 1=record: ',scatvars.server_state
print,'chirp bandwidth (MHz): ',scatvars.chirp_bandwidth_hz/1e06
print,'antenna beamwidth, (deg): ',private_config.antenna_beam_width
print,'corner reflector radar cross section (square meters): ',private_config.corner_reflector_sigma
print,'pedestal_max_speed (deg/s): ',private_config.pedestal_max_speed
print,'pedestal height (m): ',scatvars.pedestal_height_m
file_header_size=n_tags(file_header_t, /data_length) ; BYTES
waitkey,dum  
END                             ;readheader procedure 

PRO readraw,configvars,raw,scan_index,sweep_count,transition_flag,elevation,n_blocks,n_pol,scatvars,$
  filename,file_header_size,elapsed_time,time_sec,gps_latitude,gps_longitude,along_track_tilt,$
     cross_track_tilt,independent_sample_index,distance,az_proc_index,sweep_count_override



;determine number of blocks in file

  @meta_header.txt    

;@meta_header.txt

meta_header_size=n_tags(meta_header_t, /data_length) ; BYTES
print,'meta header size: ',meta_header_size
n_bytes_per_samp=2

n_pol=6;6 data products per block

block_size=scatvars.n_gates*n_pol*n_bytes_per_samp;

filedat=file_info(filename)

file_size=filedat.size
;number of data blocks:
n_blocks=(file_size-file_header_size)/(meta_header_size+block_size)
raw = fltarr(n_blocks,n_pol,scatvars.n_gates)  
time_sec=dblarr(n_blocks)
twelve_volts=fltarr(n_blocks)   
input_twelve_volts=fltarr(n_blocks)   
five_volts=fltarr(n_blocks)
minus_five_volts=fltarr(n_blocks)
RF_unit_plate_temp=dblarr(n_blocks)
LCD_display_temp=dblarr(n_blocks)
power_supply_plate_temp=dblarr(n_blocks)
power_supply_plate_temp=dblarr(n_blocks)
elevation = dblarr(n_blocks)
azimuth = dblarr(n_blocks)
sweep_count=intarr(n_blocks)
scan_index=intarr(n_blocks)
transition_flag=bytarr(n_blocks)
along_track_tilt=fltarr(n_blocks)
cross_track_tilt=fltarr(n_blocks)
block_to_block_delta=fltarr(n_blocks,6)
gps_latitude=dblarr(n_blocks)
gps_longitude=dblarr(n_blocks)
plo_20_ghz_lock=bytarr(n_blocks)

 

if configvars.i_raw_plot ge 1 then window,0,re=2
if configvars.i_raw_plot ge 1 then wshow,0
if configvars.i_corner_process eq 1 then begin
      print,'i_corner_process=1 assuming you are processing a corner reflector file'
      print,'change i_corner_process to 0 in Ku/Ka_SCAT_IDL_config.txt to process regular data'
      waitkey,dum
endif

   FOR i=0,n_blocks-1 DO BEGIN 
      readu,5,meta_header_t   
      time_sec(i)=double(meta_header_t.timestamp_seconds)+double(meta_header_t.timestamp_nanoseconds)/1e09      
      if i eq 0 then begin
         print,'start time: ',meta_header_t.timestamp_seconds, ' seconds'
         print,'nano-seconds: ',meta_header_t.timestamp_nanoseconds
      endif

      RF_unit_plate_temp(i)=meta_header_t.project_status.rcb_status.rf_plate_temp
      LCD_display_temp(i)=meta_header_t.project_status.rcb_status.lcd_display_temp
      power_supply_plate_temp(i)=meta_header_t.project_status.rcb_status.power_supply_plate_temp
      azimuth(i) = meta_header_t.project_status.ped_status.az_pos
      elevation(i) = meta_header_t.project_status.ped_status.el_pos
      sweep_count(i)=meta_header_t.project_status.ped_status.sweep_count
      scan_index(i)=meta_header_t.project_status.scan_index
      transition_flag(i)=meta_header_t.project_status.ped_status.transition_flag
      along_track_tilt(i)=meta_header_t.project_status.rcb_status.along_track_tilt
      cross_track_tilt(i)=meta_header_t.project_status.rcb_status.cross_track_tilt
      gps_longitude(i)=meta_header_t.project_status.gps_longitude
      gps_latitude(i)=meta_header_t.project_status.gps_latitude 
    plo_20_ghz_lock(i)=meta_header_t.project_status.rcb_status.plo_20_ghz_lock;only present in Ka-band system
      five_volts(i)=meta_header_t.project_status.rcb_status.five_volts
      minus_five_volts(i)=meta_header_t.project_status.rcb_status.minus_five_volts
      twelve_volts(i)=meta_header_t.project_status.rcb_status.twelve_volts
      input_twelve_volts(i)=meta_header_t.project_status.rcb_status.input_twelve_volts

      if (i/100)*100 eq i then print,'reading block number: ',i,' time: ',cmsystime(time_sec(i))


;order of data gathering: HH, VV, HV, VH, cal, noise

      hh = intarr(block_size/2/n_pol)
      vv = intarr(block_size/2/n_pol)
      hv = intarr(block_size/2/n_pol)
      vh = intarr(block_size/2/n_pol)
      cal = intarr(block_size/2/n_pol)
      noise = intarr(block_size/2/n_pol)

      readu,5,hh,vv,hv,vh,cal,noise
      
;      cal=cal+sin(findgen(block_size/2/n_pol)*2.5)*stddev(cal)

      raw(i,0,*) = counts_to_voltage(scatvars,vv)
      raw(i,1,*) = counts_to_voltage(scatvars,hv)
      raw(i,2,*) = counts_to_voltage(scatvars,vh)
      raw(i,3,*) = counts_to_voltage(scatvars,hh)
      raw(i,4,*) = counts_to_voltage(scatvars,cal)
      raw(i,5,*) = counts_to_voltage(scatvars,noise)

      pol_label=strarr(6)
      pol_label(0)='raw data, VV'
      pol_label(1)='raw data, HV'
      pol_label(2)='raw data, VH'
      pol_label(3)='raw data, HH'
      pol_label(4)='raw data, CAL'
      pol_label(5)='raw data, noise'
;if configvars.i_raw_plot ge 1 and (i/3)*3 eq i then begin
if configvars.i_raw_plot ge 1 and i eq 0 then begin
   ymax=max(raw)
   ymin=min(raw)
print,'block: ',i
   for jj=0,n_pol-1 do begin 
 ;  for jj=3,3 do begin 
               plot,raw(i,jj,*),ys=1,yr=[ymin,ymax],title=pol_label(jj),charsize=1.5,xtitle='sample number',ytitle='voltage'
               oplot,raw(i,jj,*),color=160
            print,pol_label(jj)
            waitkey,dum
         endfor
      endif
endfor;i 0 to n_blocks-1                          ;

index=where(time_sec ne 0,count)
index0=index(0:count-2)
elapsed_time=time_sec(index)-time_sec(0)
if configvars.i_temp_time_plot eq 1 then begin
window,0,retain=2,xsize=550,ysize=300,xpos=140,ypos=125
window,1,retain=2,xsize=550,ysize=300,xpos=700,ypos=125
window,2,retain=2,xsize=550,ysize=300,xpos=140,ypos=450
window,3,retain=2,xsize=550,ysize=300,xpos=700,ypos=450

wset,0
plot,elapsed_time/60,twelve_volts,title='voltages',charsize=1.5,xtitle='elapsed time (m)',ytitle='volts',xs=1,yr=[-15,15],ys=1
oplot,elapsed_time/60,twelve_volts,color=160
oplot,elapsed_time/60,input_twelve_volts,line=2,color=160
oplot,elapsed_time/60,five_volts,color=80
oplot,elapsed_time/60,-minus_five_volts,color=40
wset,1
plot,elapsed_time/60.,RF_unit_plate_temp(index),charsize=1.5,ytitle='temperature (deg. C)',title='RF unit plate temperature',xtitle='time (m)',xs=1
wset,2
plot,elapsed_time/60.,power_supply_plate_temp(index),charsize=1.5,ytitle='temperature (deg. C)',title='power supply plate temperature',xtitle='time (m)',xs=1
wset,3
plot,elapsed_time/60.,LCD_display_temp(index),charsize=1.5,ytitle='temperature (deg. C)',title='LCD display temperature',xtitle='time (m)',xs=1
waitkey,dum
 wset,0
plot,elapsed_time/60.,cross_track_tilt,charsize=1.5,ytitle='angle (deg)',title='cross-track tilt (positioner el angle ---)',xtitle='time (m)',xs=1
oplot,elapsed_time/60,elevation,line=2,color=160
wset,1
plot,elapsed_time/60.,along_track_tilt,charsize=1.5,ytitle='angle (deg)',title='along-track tilt',xtitle='time (m)',xs=1
wset,2
ymean=(long(100*mean(gps_longitude))/100.)
ymin=ymean-.05
ymax=ymean+.05
plot,elapsed_time/60.,gps_longitude,charsize=1.5,ytitle='angle (deg)',title='GPS longitude',xtitle='time (m)',xs=1,yr=[ymin,ymax],ys=1
wset,3
ymean=(long(100*mean(gps_latitude))/100.)
ymin=ymean-.05
ymax=ymean+.05
plot,elapsed_time/60.,gps_latitude,charsize=1.5,ytitle='angle (deg)',title='GPS latitude',xtitle='time (m)',xs=1,yr=[ymin,ymax],ys=1
print
print,'time between records: ',elapsed_time(1)-elapsed_time(0)
waitkey,dum
endif
close,5
independent_sample_index=indgen(n_blocks)
if max(scan_index) ne 1 then compute_distance,gps_latitude,gps_longitude,elapsed_time,scatvars,configvars,distance,independent_sample_index
if configvars.i_az_override eq 1 then select_index_az,configvars,azimuth,az_proc_index,sweep_count,sweep_count_override,elapsed_time

end;                             ;readraw procedure


PRO compute_range_profiles,range_gate_spacing,n_blocks,n_pol,ngates,gate_offset,gate_plot_max,elevation,scan_index,$
sweep_count,raw,pos_height,height,range,spec,gate_max_cal,geometric_peak,scatvars,$
private_config,configvars
if configvars.i_raw_plot ge 1 then begin
window,0,retain=2,xsize=550,ysize=300,xpos=140,ypos=125
window,1,retain=2,xsize=550,ysize=300,xpos=700,ypos=125
;window,2,retain=2,xsize=850,ysize=500,xpos=140,ypos=450
window,2,retain=2,xsize=550,ysize=300,xpos=140,ypos=450
window,3,retain=2,xsize=550,ysize=300,xpos=700,ypos=450
endif

max_display_range=configvars.max_display_range
;define constants

   pi = !pi


;compute complex range profiles
   
   raw_length=n_elements(raw(0,0,*))
   fac=alog10(raw_length)/alog10(2)
   print,'raw length',raw_length
   print,'fac: ',fac
                                ;the following code is also
                                ;implemented by the server when
                                ;finding the total corner reflector
                                ;power with ipad=1
   upfac=2
   if raw_length gt 4096  then upfac=1
   ipad=1;set to 1 to use zero padding
   if ipad eq 1 then begin
      fftlen=long(2.^(floor(fac)+upfac)) 
      print,'fft length: ',fftlen
      sizepad=fftlen-raw_length
      zeros=fltarr(sizepad)
      padfac=float(fftlen)/raw_length ;used for digital rx block-averaging equalization 
   endif else begin
      fftlen=raw_length
      padfac=1
      print,'NOTICE: no zero padding.  Just use for testing'
      print,'enter .c to continue'
      stop
   endelse



if configvars.instrument_name eq 'Ku-Scat' then begin
   range_offset_hh=-.01
   range_offset_vh=0
   range_offset_hv=0
   ped_el_offset_inches=12.6065
   hyp_inches=46
   hyp_phi_deg=0
endif
if configvars.instrument_name eq 'Ka-Scat' then begin
   range_offset_hh=0
   range_offset_vh=0
   range_offset_hv=0

   ped_el_offset_inches=12.6065
   hyp_inches=46
   hyp_phi_deg=0

endif

   range_gate_spacing=scatvars.range_resolution_m*raw_length/fftlen
   radar_range_offset=private_config.range_to_antenna_m
   gate_offset = fix(radar_range_offset/range_gate_spacing)


                                ;important to line up signals at all polarizations
                                ;to get proper correlation
                                ;coefficient   
;
;
;read these in from configuration file 
;
   
   gate_offset_hv=round(range_offset_hv/range_gate_spacing)
   gate_offset_vh=round(range_offset_vh/range_gate_spacing)
   gate_offset_hh=round(range_offset_hh/range_gate_spacing)   

   print,'gate_offset_hv: ',gate_offset_hv
   print,'gate_offset_vh: ',gate_offset_vh
   print,'gate_offset_hh: ',gate_offset_hh

   ngates=fftlen/2
   gate_plot_max = ngates-1
  
  print,'raw range resolution, range_gate_spacing with zero pad: ',scatvars.range_resolution_m,range_gate_spacing
  spec =dcomplexarr(n_blocks,n_pol,fftlen)

   range = findgen(ngates)*range_gate_spacing-radar_range_offset
   print,'maximum range (m): ',max(range)
   if max_display_range gt max(range) then max_display_range=max(range)
      itest=where(max_display_range lt range,countr)
      if countr gt 0 then gate_plot_max=itest(0)
print,'maximum displayed range (m): ',max_display_range
  weights = hanning(scatvars.n_gates)
  ; weights(*)=1
  if configvars.i_sub_band eq 1 then begin
     sub_band_fac=configvars.sub_bandwidth/(scatvars.chirp_bandwidth_hz/1e06)
     weights_sub=hanning(scatvars.n_gates*sub_band_fac)
     sweep_start_freq=(scatvars.chirp_center_frequency_hz-scatvars.chirp_bandwidth_hz/2)/1e6
     sweep_end_freq=sweep_start_freq+scatvars.chirp_bandwidth_hz/1e06
     sweep_bandwidth=scatvars.chirp_bandwidth_hz/1e06
     sub_band_start_freq=configvars.sub_center_frequency-configvars.sub_bandwidth/2.
     weights_start_index=fix(scatvars.n_gates*(sub_band_start_freq-sweep_start_freq)/(sweep_bandwidth))
     weights_end_index=fix(weights_start_index+n_elements(weights_sub))-1
     weights=fltarr(scatvars.n_gates)
     weights(weights_start_index:weights_end_index)=weights_sub
  endif

   range_index = intarr(n_blocks)

;set up variables to compute antenna height

   ped_el_offset_m=ped_el_offset_inches*.0254  
   hyp_m=hyp_inches*.0254
   hyp_phi=hyp_phi_deg*!dtor


spec_coh_ave =complexarr(n_blocks,n_pol,fftlen)
spec_coh_ave0=complexarr(n_pol,fftlen)

print,'number of records: ',n_blocks

if scatvars.decimation_mode eq 0 then equalize,scatvars,padfac,H_k_norm else H_k_norm=1.0

;compute antenna phase center height
ped_elevation_deg=elevation
pos_height=scatvars.pedestal_height_m
height=pos_height+ped_el_offset_m+hyp_m*sin(ped_elevation_deg*!dtor-hyp_phi)
geometric_peak = height/cos(elevation*!dtor) ;peak range determined from geometry

   IF gate_plot_max GT ngates-1 THEN gate_plot_max = ngates-1
i_suppress=0
sweep_count_old=0
print,'processing data into range profiles: '
if configvars.i_raw_plot ne 0 then begin
for i=0,3 do begin
wshow,i
endfor
endif
   

z_ohms=50.0                     ;
                      ;pscale preserves peak power
Watts_to_mW=1000.
pscale=Watts_to_mW*2.0*(fftlen/total(weights))^2/z_ohms
vscale=sqrt(pscale)


   FOR i=0,n_blocks-1 DO BEGIN
      if (i/100)*100 eq i then print,'processing record number: ',i
      FOR j=0,3 DO BEGIN
         if ipad eq 1 then rawpad=[weights*reform(raw(i,j,*)),zeros] else rawpad=weights*reform(raw(i,j,*))
         spec(i,j,*) =(fft(rawpad,-1)-spec_coh_ave(i,j,*))*vscale/H_k_norm ;subtract sky noise from ground data and equalize result to account for block averaging in digital rx
     ENDFOR 


      ;correct for fine offset in range 
      spec(i,1,*)=shift(spec(i,1,*),gate_offset_hv)
      spec(i,2,*)=shift(spec(i,2,*),gate_offset_vh)
      spec(i,3,*)=shift(spec(i,3,*),gate_offset_hh)
      FOR j=4,5 DO BEGIN;don't do coherent subtraction on cal signal and noise
         if ipad eq 1 then rawpad=[weights*reform(raw(i,j,*)),zeros] else rawpad=weights*reform(raw(i,j,*))
         spec(i,j,*) =fft(rawpad,-1)*vscale/H_k_norm
      ENDFOR;
      IF configvars.i_raw_plot EQ 1 and scan_index(i) eq 1 THEN BEGIN 
         IF sweep_count(i) NE sweep_count_old THEN BEGIN 
            plot_range_profiles,i,spec,sweep_count,ngates,gate_offset,gate_plot_max,geometric_peak,elevation,range
         ENDIF 
      ENDIF
      if configvars.i_raw_plot eq 2 and abs(scan_index(i)) eq 1 or configvars.i_raw_plot eq 3 then plot_range_profiles,i,spec,sweep_count,ngates,gate_offset,gate_plot_max,geometric_peak,elevation,range
      sweep_count_old=sweep_count(i)
   ENDFOR;i loop
   
   pspec_cal = 20*alog10(abs(spec(0,4,gate_offset:ngates-1)))
   gate_max_cal = gate_offset+where(pspec_cal EQ max(pspec_cal))
   index_surface_vec=where(scan_index eq 1,count)
   index_surface=index_surface_vec(0)
END                             ;compute_range_profiles procedure

PRO plot_range_profiles,i,spec,sweep_count,ngates,gate_offset,gate_plot_max,geometric_peak,elevation,range
;

   ymin = -90
   ymax = 10
   yoff=-20
   yoff2=-25
   pspec_vv = 20*alog10(abs(spec(i,0,0:ngates-1)))
   pspec_hh = 20*alog10(abs(spec(i,3,0:ngates-1)))
;check on ordering of subscripts 
   pspec_vh = 20*alog10(abs(spec(i,1,0:ngates-1)))
   pspec_hv  = 20*alog10(abs(spec(i,2,0:ngates-1)))
   pspec_cal = 20*alog10(abs(spec(i,4,0:ngates-1)))
   pspec_noise = 20*alog10(abs(spec(i,5,0:ngates-1)))
   phase_vv=atan(imaginary(spec(i,0,0:ngates-1)),float(spec(i,0,0:ngates-1)))/!dtor
   phase_hv=atan(imaginary(spec(i,1,0:ngates-1)),float(spec(i,1,0:ngates-1)))/!dtor
   phase_vh=atan(imaginary(spec(i,2,0:ngates-1)),float(spec(i,2,0:ngates-1)))/!dtor
   phase_hh=atan(imaginary(spec(i,3,0:ngates-1)),float(spec(i,3,0:ngates-1)))/!dtor

   phase_cal=atan(imaginary(spec(i,4,0:ngates-1)),float(spec(i,4,0:ngates-1)))/!dtor
   phase_noise=atan(imaginary(spec(i,5,0:ngates-1)),float(spec(i,5,0:ngates-1)))/!dtor
   phase_vv_hh=phase_vv-phase_hh

   indexvvhh=where(phase_vv_hh lt -180.0,countvvhh)
   if countvvhh gt 0 then phase_vv_hh(indexvvhh)=phase_vv_hh(indexvvhh)+360.0
   indexvvhh=where(phase_vv_hh gt 180.0,countvvhh)
   if countvvhh gt 0 then phase_vv_hh(indexvvhh)=phase_vv_hh(indexvvhh)-360.0
;print,'***********************************'
;   print,'max pspec_vv: ',max(pspec_vv(gate_offset:gate_plot_max))
;   print,'max pspec_vh: ',max(pspec_vh(gate_offset:gate_plot_max))
;   print,'max pspec_hv: ',max(pspec_hv(gate_offset:gate_plot_max))
;   print,'max pspec_hh: ',max(pspec_hh(gate_offset:gate_plot_max))
;   print,'max pspec_cal: ',max(pspec_cal(gate_offset:gate_plot_max))
;   print,'difference, max VH max cal: ',max(pspec_vh(gate_offset:gate_plot_max))-max(pspec_cal(gate_offset:gate_plot_max))
;print,'***********************************'
   
;printf,1,max(pspec_vh(gate_offset:gate_plot_max))
;printf,1,max(pspec_cal(gate_offset:gate_plot_max))
;indexVH=gate_offset+where(pspec_vh eq max(pspec_vh(gate_offset:gate_plot_max)))
;indexcal=gate_offset+where(pspec_cal eq max(pspec_cal(gate_offset:gate_plot_max)))
;print,'***********************************'
;pspec_vh_sum=10*alog10(total(10^(pspec_vh(indexvh-4:indexvh+4)/10.)))
;pspec_cal_sum=10*alog10(total(10^(pspec_cal(indexcal-4:indexcal+4)/10.)))
;   print,'total pspec_vh: ',pspec_vh_sum
;   print,'total pspec_cal: ',pspec_cal_sum
;   print,'difference: ',pspec_vh_sum-pspec_cal_sum
;print,'***********************************'

   wset,2
   wshow,2
   plot,range(gate_offset:gate_plot_max),pspec_vv(gate_offset:gate_plot_max),xtitle='range (m)',ytitle='power (dB)',title='range profiles VV',charsize=1.6,yrange=[ymin,ymax],xs=1,ys=1
   oplot,range(gate_offset:gate_plot_max),pspec_vv(gate_offset:gate_plot_max),color=160
   wset,3
   wshow,3
   plot,range(gate_offset:gate_plot_max),pspec_hh(gate_offset:gate_plot_max),xtitle='range (m)',ytitle='power (dB)',title='range profiles HH',charsize=1.6,yrange=[ymin,ymax],xs=1,ys=1
  oplot,range(gate_offset:gate_plot_max),pspec_hh(gate_offset:gate_plot_max),color=80
  wset,0
  wshow,0
   plot,range(gate_offset:gate_plot_max),pspec_hv(gate_offset:gate_plot_max),xtitle='range (m)',ytitle='power (dB)',title='range profiles HV/VH',charsize=1.6,yrange=[ymin,ymax],xs=1,ys=1
  oplot,range(gate_offset:gate_plot_max),pspec_hv(gate_offset:gate_plot_max),color=60
  oplot,range(gate_offset:gate_plot_max),pspec_vh(gate_offset:gate_plot_max),color=40
  wset,1
  wshow,1
   plot,range(gate_offset:gate_plot_max),pspec_vv(gate_offset:gate_plot_max),xtitle='range (m)',ytitle='power (dB)',title='range profiles VV/HH/VH',charsize=1.6,yrange=[ymin,ymax],xs=1,ys=1
   oplot,range(gate_offset:gate_plot_max),pspec_vv(gate_offset:gate_plot_max),color=160
  oplot,range(gate_offset:gate_plot_max),pspec_hh(gate_offset:gate_plot_max),color=80
  oplot,range(gate_offset:gate_plot_max),pspec_vh(gate_offset:gate_plot_max),color=40
  if elevation(i) gt 0.0 and elevation(i) le 89.50 then begin
     rvec = [geometric_peak(i),geometric_peak(i)]
     yvec = [ymin,ymax]
     oplot,rvec,yvec,color=1000
   xyouts,.85*range(gate_plot_max),yoff2,'elevation angle (deg): ',alignment=1,charsize=1.4
   xyouts,.85*range(gate_plot_max),yoff2,strtrim(string(round(elevation(i))),2),charsize=1.4
   endif
;cal signal
waitkey,dum
;goto,skipcalplot
wset,0
wshow,0
plot,range(gate_offset:gate_plot_max),pspec_cal(gate_offset:gate_plot_max),charsize=1.5,xtitle='range(m)',ytitle='power (dB)',title='internal cal loop',yr=[ymin,ymax],ys=1,xs=1
  oplot,range(gate_offset:gate_plot_max),pspec_cal(gate_offset:gate_plot_max),color=160

goto,skip_vv_zoom
wset,1
wshow,1
ymax=0
ymin=-12
plot,range(gate_offset:gate_plot_max),pspec_vv(gate_offset:gate_plot_max)-max(pspec_vv(gate_offset:gate_plot_max)),charsize=1.5,xtitle='range(m)',ytitle='power (dB)',title='VV power/internal cal loop',yr=[ymin,ymax],ys=1,xs=1,xr=[1.5,2.0]
  oplot,range(gate_offset:gate_plot_max),pspec_vv(gate_offset:gate_plot_max)-max(pspec_vv(gate_offset:gate_plot_max)),color=160
  oplot,range(5*gate_offset:gate_plot_max),pspec_cal(5*gate_offset:gate_plot_max)-max(pspec_cal(gate_offset:gate_plot_max)),color=80
skip_vv_zoom:  
waitkey,dum
skipcalplot:
END                             ;plot_range_profiles procedure



PRO calibrate,calvars,configvars,gate_offset,gate_plot_max,spec,ngates,n_blocks,n_pol,range,gate_max_cal,reference_calibration_loop_power,$
current_calibration_loop_power,ainv,finv,range_gate_spacing,corner_range,scatvars,scan_index,sweep_count,$
elevation,total_corner_power_vv,total_corner_power_hh,base_filename,calfile

;
;find internal calibration loop power
;
  calpwr_i =abs(spec(0:n_blocks-1,4,gate_max_cal))^2 
  median_cal_pwr = median(calpwr_i)
  calpwr_i_db = 10*alog10(calpwr_i)
  median_cal_pwr_db = 10*alog10(median_cal_pwr)
  current_calibration_loop_power = total(calpwr_i)/n_blocks
  print,'mean calibration power, dB',10*alog10(current_calibration_loop_power)
  if configvars.i_calibrate eq 0 then begin               ; using old cal file
    read_calfile,ainv,finv,reference_calibration_loop_power,corner_range,total_corner_power_vv,total_corner_power_hh,calfile
print,'mean calibration power stored in calibration file, dB',10*alog10(reference_calibration_loop_power)
ENDIF

if configvars.i_calibrate eq 1 then begin
;uses clutter calibration procedure to find fmat and amat.  Assumes
;neglible antenna cross-polarization and isotropic, monostatic
;scattering from surface or volume target such that angle of
;Svv*conj(Shh)=0 and angle of Svh*conj(Shv)=0.

   amat_fmat_compute,spec,scatvars,configvars,calvars,range_gate_spacing,n_blocks,range,ngates,ainv,finv,corner_range,total_corner_power_vv,total_corner_power_hh,scan_index,sweep_count,elevation,current_calibration_loop_power,base_filename,calfile
   
   reference_calibration_loop_power=current_calibration_loop_power
endif

END                             ;calibrate procedure

PRO polsig,L_matrix,iline,control_plots,line_el_bin,i_multi_bin
  wset,0
                                ;line_el_bin is either the line elevation or the range bin index
  
;this procedure plots polarization signatures which are a graphical
;representation of the response of the target to prescribed
;combination of transmit and receive polarizations. The polarization
;signatures computed in this subroutine are so-called "co-polarized"
;signatures, which displays the average power received by a radar illuminating
;the target with a variable polarization antenna as the polarization
;takes on all possible states.  
 
   nchi = 25; number of chi angles generated; chi is the polarization ellipse's ellipticity (e.g., chi=0 for near polarization; chi=+/-45 for RHCP/LHCP)
   npsi = 41;number of psi angles generated; psi is the orientation of the polarization ellipse (e.g., psi=0 is horizontal polarization; psi=90 is vertical polarization)
   chivec = fltarr(nchi)
   psivec = fltarr(npsi)
   M_matrix = fltarr(4,4)       ;Stokes scattering matrix

   p_rec = fltarr(nchi,npsi)
   p_total = fltarr(nchi,npsi)
   p_rec_xpol = fltarr(nchi,npsi)
   degree_of_polarization=fltarr(nchi,npsi)
;   Gmat =fltarr(4)
   Gmat =fltarr(1,4)
;   Gmat_xpol =fltarr(4)
   Gmat_xpol =fltarr(1,4)

      
   M_matrix=L_matrix
;   M_matrix(3,*) = -M_matrix(3,*); M_matrix is the Stokes Scattering
;   matrix, similar to the Mueller matrix, with a sign change on the
;   last row. 

   M_matrix(*,3) = -M_matrix(*,3); M_matrix is the Stokes Scattering matrix, similar to the Mueller matrix, with a sign change on the last row. 
      
      FOR k=0,nchi-1 DO BEGIN 
         FOR j=0,npsi-1 DO BEGIN 
            chideg = -45.+90*k/(nchi-1)
            psideg = -90.+180.*j/(npsi-1) ; vpol has psi=0 deg; hpol, psi=+/-90.0 degrees. See Yamaguchi and Boerner, pg 174 SPIE Vol 1748 Radar Polarimetry 1992 with vpol aligned with X axis.
            chi = chideg*!dtor
            psi = psideg*!dtor
            
            chi_xpol=-chi
            psideg_xpol=psideg+90.0
            index=where(psideg_xpol gt 90.0,count)
            if count gt 0 then psideg_xpol(index)=psideg_xpol-180.0
            psi_xpol=psideg_xpol*!dtor
           
            ;gmat is the Stokes vector describing the antenna polarization
            gmat(0,0) = 1.0
            gmat(0,1) = cos(2.*psi)*cos(2.*chi)
            gmat(0,2) = sin(2.*psi)*cos(2.*chi)
            gmat(0,3) = sin(2.*chi)

                                ;gmat_xpol is the Stokes vector of a
                                ;cross-polarized antenna relative to
                                ;the copol antenna state 
            gmat_xpol(0,0) = 1.0
            gmat_xpol(0,1) = cos(2.*psi_xpol)*cos(2.*chi_xpol)
            gmat_xpol(0,2) = sin(2.*psi_xpol)*cos(2.*chi_xpol)
            gmat_xpol(0,3) = sin(2.*chi_xpol)

            MG = M_matrix##gmat; is the Scattered Stokes Vector
            p_rec(k,j)= total(gmat*MG);gmat*MG is the dot-product of the scattered Stokes vector and the antenna polarization; this dot product yields the received power
          p_rec_xpol(k,j)= total(gmat_xpol*MG);gmat*MG is the dot-product of the scattered Stokes vector and the antenna polarization; this dot product yields the received power
          degree_of_polarization(k,j)=sqrt(MG(1)^2+MG(2)^2+MG(3)^2)/MG(0)
          
          p_total(k,j)=2*MG(0)
             chivec(k) = chideg
            psivec(j) = psideg
            
         ENDFOR 
      ENDFOR
      print
      if i_multi_bin eq 0 and line_el_bin ne 0 then  print,'normalized Mueller matrix for scan line: ', ' = ',2*iline 
      if i_multi_bin eq 0 and line_el_bin eq 0 then  print,'normalized Mueller matrix averaged over all bins: ' 
      if i_multi_bin eq 1 then  print,'normalized Mueller matrix for bin number: ', ' = ',line_el_bin 
      print,L_matrix/max(L_matrix)
      wset,0
      plot_ptotal = p_total/max(p_total)
      plot_rec = p_rec/max(p_total)
      plot_rec_xpol = p_rec_xpol/max(p_total)
;print,'max plot_rec,plot_rec_xpol: ',max(plot_rec),max(plot_rec_xpol)
;if max(plot_rec) gt 1 or max(plot_rec_xpol) gt 1 then stop

L_norm=L_matrix/max(L_matrix)
LDR_i=(1-L_norm(1,1))/(1+L_norm(1,1))
print,'LDR (dB): ',10*alog10(LDR_i)

;      print,'max normalized cross-pol response: ',max(plot_rec_xpol)
if i_multi_bin eq 0 then titlestr='co-pol el='+strmid(strtrim(string(line_el_bin),2),0,3) else titlestr='co-pol range bin='+strmid(strtrim(string(line_el_bin),2),0,3)
      surface,plot_rec,chivec,psivec,ax=40,az=40.,title=titlestr,$
ztitle='relative NRCS',ytitle='psi',xtitle='chi',charsize=3,yticks=3,ytickv=[-90,0,90],xticks=3,xtickv=[-45.,0,45],zrange=[0,1]
if control_plots eq 1 then begin
   wset,1
   if i_multi_bin eq 0 then titlestr='cross-pol el='+strmid(strtrim(string(line_el_bin),2),0,3) else titlestr='cross-pol range bin='+strmid(strtrim(string(line_el_bin),2),0,3)
      surface,plot_rec_xpol,chivec,psivec,ax=40,az=40.,title=titlestr,ztitle='relative NRCS',ytitle='psi',xtitle='chi',charsize=3,yticks=3,ytickv=[-90,0,90],xticks=3,xtickv=[-45.,0,45],zrange=[0,1]
      wset,2
if i_multi_bin eq 0 then titlestr='deg. of polar. el='+strmid(strtrim(string(line_el_bin),2),0,3) else titlestr='deg. of polar. range bin='+strmid(strtrim(string(line_el_bin),2),0,3)
      surface,degree_of_polarization,chivec,psivec,ax=40,az=40.,$
title=titlestr,ztitle='deg. of pol.',ytitle='psi',xtitle='chi',charsize=3,yticks=3,ytickv=[-90,0,90],xticks=3,xtickv=[-45.,0,45],zrange=[0,1]

wset,3
if i_multi_bin eq 0 then titlestr='total power el='+strmid(strtrim(string(line_el_bin),2),0,3) else titlestr='total power range bin='+strmid(strtrim(string(line_el_bin),2),0,3)
      surface,plot_ptotal,chivec,psivec,ax=40,az=40.,$
title=titlestr,ztitle='deg. of pol.',ytitle='psi',xtitle='chi',charsize=3,yticks=3,ytickv=[-90,0,90],xticks=3,xtickv=[-45.,0,45],zrange=[0,1]
endif
print,'max degree of polarization: ',max(degree_of_polarization)
print,'max copol power: ',max(plot_rec)
print,'min copol power: ',min(plot_rec)
print,'max xpol power: ',max(plot_rec_xpol)
print,'min xpol power: ',min(plot_rec_xpol)

waitkey,dum

END
PRO polscatter,smat,npoints,line_el
;this procedure provides a 2-D histogram of the polarization state of the instantaneous response to vertical and horizontal transmit polarization. 
;Change evec to find response to other polarization pairs
;This is useful for visualizing how the target depolarized the
;incident wave; A highly polarized target, such as the corner
;reflector has a very narrow distribution of polarization states
;whereas a highly depolarizing target shows a wide distribution
;
   smat0 = reform(smat,2,2,npoints)
   evec = complexarr(2,2)

;first transmit polarization state (vertical polarization is evec(0,0)=1.0; evec(1,0).0); for lhcp: evec(0,0)=1.0; evec(1,0)= complex(0.0,1.0) 
   evec(0,0) = 1.0
   evec(1,0) = 0.0              ; 

;second transmit polarization state (horizontal polarization has evec(0,1)=0.0; evec(1,1)=1.0)
   evec(0,1) = 0.0
   evec(1,1) = 1.0

   escat = complexarr(2,2,npoints) ;vector of scattered polarization states

;   escat format:  [Ev1   Eh1]
;                  |         |
;                  [Ev2   Eh2]
;
;   where first subscript is the rx polarization state, second: tx polarization number (first or second tx polarization state)
;    

   FOR k=0,npoints-1 DO BEGIN 
      escat(*,*,k) = smat0(*,*,k)#evec
   ENDFOR 

   FOR i=0,1 DO BEGIN ; i=0 for first transmit polarization; i=1 for second transmit polarization
      magv = abs(escat(0,i,*))
      magh = abs(escat(1,i,*))
      alpha_h = atan(imaginary(escat(1,i,*)),real_part(escat(1,i,*)))
      alpha_v = atan(imaginary(escat(0,i,*)),real_part(escat(0,i,*)))
      delta = alpha_h-alpha_v ;phase difference between v and h field components
      psi= .5*atan((2*magv*magh*cos(delta)),(magv^2-magh^2)) ;polarization ellipse orientation angle
      chi = .5*asin((2*magv*magh*sin(delta))/(magv^2+magh^2)) ;polarization ellipse ellipticity angle
      density = hist_2d(chi*57.3,psi*57.3,bin1=3,bin2=3,min1=-45,max1=45,min2=-90,max2=90) ;histogram of psi,chi
      szd = size(density)
      nchi = szd(1)
      npsi = szd(2)
;angle vectors for labelling axes 
      chivec = -45.+90*findgen(nchi)/(nchi-1)
      psivec = -90.+180.*findgen(npsi)/(npsi-1)
   
   IF i EQ 0 THEN begin
       wset,2
       wshow,2
       surface,density,chivec,psivec,ax=40,az=40.,title='V tx histogram el='+$
strmid(strtrim(string(line_el),2),0,2),ztitle='number of hits',ytitle='psi, deg',$
xtitle='chi, deg',charsize=3,yticks=3,ytickv=[-90,0,90],xticks=3,xtickv=[-45.,0,45]
   endif

      IF i EQ 1 THEN begin
          wset,3
          wshow,3
          surface,density,chivec,psivec,ax=40,az=40.,title='H tx histogram el='+strmid(strtrim(string(line_el),2),0,2),ztitle='number of hits',ytitle='psi, deg',xtitle='chi, deg',charsize=3,yticks=3,ytickv=[-90,0,90],xticks=3,xtickv=[-45.,0,45]
      endif

;      print,'enter "s" to store .eps of histogram files'
      dum = 0
;      waitkey,dum
      GOTO,jumpplot
      IF dum EQ 's' THEN BEGIN 
         filepol = strarr(1)
         IF i EQ 0 THEN read,'enter file name for vertical tx (.eps): ',filepol
         IF i EQ 1 THEN read,'enter file name for horizontal tx (.eps): ',filepol
         set_plot,'ps'
         device,file=filepol,/helvetica,yoffset=8,/encapsulated
         IF i EQ 0 THEN surface,density,chivec,psivec,ax=30,az=35.,ztitle='number of hits',ytitle='psi, deg',xtitle='chi, deg',charsize=3,yticks=5,ytickv=[-90,-45,0,45,90],xticks=3,xtickv=[-45.,0,45],xthick=4,ythick=4,charthick=4
         IF i EQ 1 THEN surface,density,chivec,psivec,ax=30,az=35.,ztitle='number of hits',ytitle='psi, deg',xtitle='chi, deg',charsize=3,yticks=5,ytickv=[-90,-45,0,45,90],xticks=3,xtickv=[-45.,0,45],xthick=4,ythick=4,charthick=4
         device,/close
         set_plot,'x'
      ENDIF 
jumpplot:
   ENDFOR 
waitkey,dum
END;pol_scatter


PRO read_calfile,ainv,finv,reference_calibration_loop_power,corner_range,total_corner_power_vv,total_corner_power_hh,calfile
  if calfile ne "select" then begin
     print,'******************'
     print,'calibration file in use: ',calfile
     print,'edit name of calibration file in: Ku-(or Ka-)SCAT_IDL_config.txt if necessary'
    ; waitkey,dum
  endif

  if calfile eq "select" then calfile=pickfile(file="",path=".",/read,filter="*.cal",title="Pick a file to open")
     close,1
     openr,1,calfile
     ainv = complexarr(2,2)
     finv = complexarr(2,2)
     reference_calibration_loop_power=fltarr(1)
     corner_range=0.0
     total_corner_power_vv=0.0
     total_corner_power_hh=0.0
     amat=complexarr(2,2)
     fmat=complexarr(2,2)
     readf,1,ainv,finv,reference_calibration_loop_power,corner_range,total_corner_power_vv,total_corner_power_hh,amat,fmat
     close,1
    ; amat = invert(ainv)
    ; fmat = invert(finv)
     vmat_cr = amat#fmat
     print
     print,'receiver distortion matrix from cal file: '
     print,amat
     print
     print,'cal file transmitter distortion matrix from cal file: '
     print,fmat
     print
     print,'reference cal loop power: ',reference_calibration_loop_power
     print
     print,'corner range: ',corner_range
     print
     print,'integrated corner reflector power, v: ',total_corner_power_vv
     print
     print,'integrated corner reflector power, h: ',total_corner_power_hh

  end;read_calfile

PRO NEXRAD_ct

;*** NEXRAD Color Table with white removed except at min value (black
;at max value)***

dc=240
h=254
;note: r,g,b values scale from lowest to highest intensity the colors
;red, green and blue
;white: r,g,b=255; black: r,g,b=0

 r=transpose([20,  20,   20,   20,   20,   20, 255, 231, 255, 255, 214, 192, 255, 153, 20])
 g=transpose([236, 160, 120, 255, 200, 144, 255, 192, 144,   20,   20,   20,   20,  85, 20])
 b=transpose([236, 246, 190,   20,   20,   20,   20,   20,   20,   20,   20,   20, 255, 201, 20])
 tr=transpose(reform([r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r],1, 240))
 tg=transpose(reform([g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g],1, 240))
 tb=transpose(reform([b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b],1, 240))

 r=bytscl(findgen(256))
 g=r
 b=r

 r(1:h)=interpolate(tr,dc*findgen(h)/(h-1))
 g(1:h)=interpolate(tg,dc*findgen(h)/(h-1))
 b(1:h)=interpolate(tb,dc*findgen(h)/(h-1))

 r(0)=255 & g(0)=255 & b(0)=255
 r(255)=0 & g(255)=0 & b(255)=0
 tvlct, r,g,b

END


PRO get_variable,filename,variable,test_string
space=''
title=''
readf,1,space
readf,1,title

string_size=strlen(strtrim(test_string))

if strmid(title,0,string_size) ne test_string then begin
print
print,'bad input from ',filename,'; expecting the following variable header: ',test_string
print
print,'instead, program read: ',strmid(title,0,string_size)
print
;stop
endif
readf,1,variable
end;get_variable

;beginning of function cmsystime
;+
; NAME:
;   CMSYSTIME
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Compute seconds since Jan 1, 1970 and (Modified) Julian Days
;
; CALLING SEQUENCE:
;   TIMEVAL1 = CMSYSTIME(TIMEVAL0, ...)
;
; DESCRIPTION: 
;
;   CMSYSTIME serves two functions.  It computes the current time in a
;   fashion similar to the built-in IDL system function SYSTIME().  It
;   also can convert between various time representations and systems,
;   including a textual format.
;
;   The current time can be obtained by invoking CMSYSTIME with the
;   /NOW keyword (which is entirely equivalent to SYSTIME(1)).
;
;   The most substantial part of CMSYSTIME, which distinguishes it
;   from SYSTIME, is its ability to convert between different time
;   formats.  CMSYSTIME recognizes can recognize and convert between
;   time in seconds (seconds since Jan 1, 1970 [ = SEC ]) and days
;   (Julian days [ = JDAY ] or "Modified" Julian days [ = MJD = JDAY -
;   2400000.5 ]).  It can also recognize and convert between local and
;   GM time.  
;
;   CMSYSTIME takes maximum care to preserve the full numerical
;   precision of the time values.  It converts all values to double
;   precision and may return days and seconds with fractional parts.
;
;   CMSYSTIME can also represent any time textually, not just the
;   current time.  The following textual formats are supported:
;        DOW MMM DD hh:mm:ss YYYY              - (Default - same as SYSTIME)
;        DOW MMM DD YYYY hh:mm:ss.uuuuuu TTTTT - (/EXTENDED)
;   where DOW and MMM are the abbreviated day of week and month in
;   English, DD is the day of the month, YYYY is the year, hh:mm:ss is
;   the time in 24 hr military time, uuuuuu are additional
;   microseconds, TTTTT is the timezone offset (in +hhmm
;   representation).
;
;   CMSYSTIME accepts one parameter, the input time to be converted.
;   Unlike SYSTIME, the *function* of CMSYSTIME is governed by various
;   keywords, as summarized in the following table:
;
;   Converting from                       Converting to
;   ---------------                       -------------
;   JDAY - /FROM_JULIAN                   JDAY - /JULIAN
;   MJD  - /FROM_MJD                      MJD  - /MJD
;   SEC  - (Default)                      SEC  - /SECONDS
;   Current time - /NOW                   TEXT - (Default or /EXTENDED)
;
;   Local time - /FROM_LOCAL              Local time - /LOCAL
;   GM time - (Default)                   GM time - (Default)
;   
;   If no argument is specified, the default is to report the current
;   time textually in the GM time zone.  CMSYSTIME automatically
;   determines the local time zone.
;
; INPUTS:
;
;   TIMEVAL0 - input time, in seconds or days, as described above.
;              This value is ignored if the NOW keyword is set.  Array
;              values are allowed.
;
; KEYWORDS:
;
;   NOW - If set, TIMEVAL0 is ignored and the current time is used as
;         input.
;
;   FROM_JULIAN - If set, TIMEVAL0 is in Julian days.
;   FROM_MJD    - If set, TIMEVAL0 is in Modified Julian days (MJD).
;   FROM_LOCAL  - If set, TIMEVAL0 is in the local time zone.
;                 If no FROM_ keywords are set, the input is assumed
;                 to be seconds from Jan 1, 1970.
;
;   JULIAN  - If set, the input is converted to Julian days upon output.
;   MJD     - If set, the input is converted to MJD upon output.
;   SECONDS - If set, the input is converted to seconds from Jan
;             1, 1970 upon output.
;   LOCAL   - If set, the input is converted to the local time zone.
;             If no "destination" keywords are set, the output is
;             converted to textual representation.
;
;   EXTENDED - Convert to a textual representation with additional
;              information, as noted above.
;
;   TIMEZONE - Upon output, the timezone offset is returned in this
;              keyword.  The offset is time difference in seconds
;              between GM time and the local time, such that LOCALTIME
;              = GMTIME + TIMEZONE
;
; RETURNS:
;   The resulting converted time(s), either as a double precision
;   number or a string.
;
; EXAMPLE:
;   
;   The equivalent to SYSTIME(0)
;     IDL> print, systime(0) & print, cmsystime(/now, /local)
;     Wed Jul  5 12:10:46 2000
;     Wed Jul  5 12:10:46 2000
;
;   The equivalent to SYSTIME(1)
;     IDL> print, systime(1) & print, cmsystime(/now,/seconds)
;        9.6277750e+08
;        9.6277750e+08
;
;   Comparison between local and GM time zones (I live in the Eastern
;    US, daylight savings)
;     IDL> print, cmsystime(/now,/extended)
;     Wed Jul  5 2000 16:13:15.659000 -0400
;     IDL> print, cmsystime(/now,/local,/extended)
;     Wed Jul  5 2000 12:13:15.664000 -0400
;    
;   What day of the week was it 200 days ago?  (Note, there are 86400
;    seconds in one day)
;     IDL> today = cmsystime(/now,/seconds)
;     IDL> print, cmsystime(today-86400L*200, /local)
;     Sat Dec 18 12:17:52 1999
;    
;
; SEE ALSO:
;
;   SYSTIME, JULDAY, CALDAT
;
; MODIFICATION HISTORY:
;   Written, CM, 05 Jul 2000
;   Printed time zone is zero when LOCAL=0, CM, 21 Aug 2000
;   Corrected behavior of /MJD (Thanks to Marshall Perrin), 03 Jun
;     2002
;   Corrected local vs. UTC problem caused by fractional UTC seconds,
;     (thanks to J. Wolfe) CM, 28 Dec 2005
;   Corrected problem with Julian day arrays, (thanks to W. Landsman),
;     CM, 29 Dec 2005
;
;  $Id: cmsystime.pro,v 1.5 2005/12/29 18:07:48 craigm Exp $
;
;-
; Copyright (C) 2000,2002,2005, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

forward_function cmsystime_xmod, cmsystime

;; Comput X MOD M, ensuring positive remainder
function cmsystime_xmod, x, m
  return, (((x MOD m) + m) MOD m)
end

;; Convert from MJD to YR/MO/DAY
pro cmsystime_mjd2ymd, mjd, yr, mo, da
  offset = 2400000.5D
  offset_int = floor(offset)         ;; Integer part of offset
  offset_fra = offset - offset_int   ;; Fractional part of offset
  nn = offset_fra + mjd
  jd_fra = cmsystime_xmod(nn+0.5D, 1D) - 0.5D
  nn = nn + offset_int - jd_fra
  nn = nn + (floor(floor((nn - 4479.5D)/36524.25D) * 0.75D + 0.5D)-37.D)
  
  yr = long(floor(nn/365.25D) - 4712.D)
  dd = floor(cmsystime_xmod(nn-59.25D, 365.25D))
  
  mo = floor(cmsystime_xmod( floor((dd+0.5D)/30.6D) + 2.D, 12.D ) + 1.D)
  da = floor(cmsystime_xmod(dd+0.5D, 30.6D) + 1.D ) + 0.5D + jd_fra
end

function cmsystime, arg0, now=now, extended=extended, $
                    local=local, from_local=from_local, $
                    julian=jul, from_julian=from_julian, $
                    mjd=mjd, from_mjd=from_mjd, $
                    seconds=seconds, timezone=timezone

  common cmsystime_common, cmsystime_timezone, cmsystime_months, $
    cmsystime_dow

  ;; Precompute names of days in week and month
  if n_elements(cmsystime_months) EQ 0 then begin
      cmsystime_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', $
                          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
      cmsystime_dow = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
  endif

  ;; Starting epoch, expressed in MJD and Julian days
  MJD_1970 = 40587D
  JD_1970  = MJD_1970 + 2400000.5D

  ;; Figure the time zone automatically, the first time around
  if n_elements(cmsystime_timezone) EQ 0 then begin
      ;; The GM time, converted to MDY
      gmtime = systime(1)
      cltime = systime(0)
      gmfrac = gmtime MOD 86400
      gm_mjd = floor(gmtime-gmfrac)/86400D + MJD_1970
      cmsystime_mjd2ymd, gm_mjd, gm_yr, gm_mo, gm_da
      gm_da = round(gm_da)

      ;; The local time
      ltime = strtrim(str_sep(strcompress(cltime),' '),2)
      da = floor(long(ltime(2)))
      ltimes = double(str_sep(ltime(3), ':'))
      lfrac = ltimes(2) + 60D*(ltimes(1) + 60D*ltimes(0))

      ;; The timezone difference...
      tz = lfrac - gmfrac
      ;; ... but we must account for day wrap-around
      if      gm_da EQ da - 1 then tz = tz + 86400 $
      else if gm_da EQ da + 1 then tz = tz - 86400 $
      else if gm_da LT da     then tz = tz - 86400 $  ;; ...and month roll-over
      else if gm_da GT da     then tz = tz + 86400    ;; ...and month roll-over

      ;; Store the new value
      cmsystime_timezone = round(tz/60)*60  ;; Round to nearest minute
  endif
  timezone = cmsystime_timezone

  ;; Compute the timezone offset, depending on which way the
  ;; conversion will go.
  offset = 0D
  if keyword_set(from_local) then offset = offset - timezone
  if keyword_set(local)      then offset = offset + timezone

  ;; Extract the time value either from the clock, or from the user
  ;; parameter
  if keyword_set(now) then begin
      ;; From clock (GMT)
      NOW_TIME:
      arg = systime(1)
      if keyword_set(from_local) then offset = 0D

  endif else begin
      ;; From user parameter
      if n_elements(arg0) EQ 0 then goto, NOW_TIME
      arg = double(arg0)

      if keyword_set(from_mjd) then begin
          ;; Convert from MJD ... avoid loss of numerical precision
          if keyword_set(mjd) then return, arg + offset
          if keyword_set(jul) then return, arg + offset + (JD_1970-MJD_1970)

          ;; Convert to seconds
          arg = (arg - MJD_1970) * 86400D

      endif else if keyword_set(from_julian) then begin
          ;; Convert from JD ... avoid loss of numerical precision if poss.
          if keyword_set(mjd) then return, arg + offset - (JD_1970-MJD_1970)
          if keyword_set(jul) then return, arg + offset

          ;; Convert to seconds
          arg = (arg - JD_1970)  * 86400D

      endif
  endelse

  ;; Add timezone offset
  if offset NE 0 then arg = arg + offset
  if keyword_set(seconds) then return, arg
  if keyword_set(jul)     then return, (arg / 86400D) +  JD_1970
  if keyword_set(mjd)     then return, (arg / 86400D) + MJD_1970

  ;; Convert to MJD, from there to MDY
  mjd = floor(arg/86400D) + MJD_1970
  dsecs = arg-floor(arg/86400D)*86400D
  hr = floor(dsecs / 3600) & dsecs = dsecs - hr*3600
  mi = floor(dsecs / 60)   & dsecs = dsecs - mi*60
  se = dsecs
  cmsystime_mjd2ymd, mjd, yr, mo, da

  ;; Day of week is simple to calculate, assumes 13 May 2000 was a Sunday
  dow = cmsystime_xmod((floor(mjd) - 51678L), 7L)

  ;; Compute the string values, unfortunately on an individual basis
  n = n_elements(yr)
  result = strarr(n)
  if keyword_set(extended) then begin
      for i = 0L, n-1 do begin
          sei = floor(se(i))
          sef = floor((se(i) - sei)*1000000D)
          result(i) = string(cmsystime_dow(dow(i)), cmsystime_months(mo(i)-1),$
                             da(i), yr(i), hr(i), mi(i), sei, sef, $
                             format=('(A3," ",A3," ",I2," ",I4.4," ",' + $
                                     'I2.2,":",I2.2,":",I2.2,".",I6.6)'))
      endfor

      ;; Extended string value includes time zone offset
      if keyword_set(local) then tzz = timezone else tzz = 0L
      tzabs = abs(tzz)
      tzhr = floor(tzabs/3600)
      tzstring = string(tzhr, floor(tzabs - tzhr*3600), $
                        format='(I2.2,I2.2)')
      if tzz LT 0 then tzstring = ' -'+tzstring $
      else             tzstring = ' +'+tzstring
      result = result + tzstring

  endif else begin
      for i = 0L, n-1 do begin
          result(i) = string(cmsystime_dow(dow(i)), cmsystime_months(mo(i)-1),$
                             da(i), hr(i), mi(i), floor(se(i)), yr(i), $
                             format=('(A3," ",A3," ",I2," ",' + $
                                     'I2.2,":",I2.2,":",I2.2," ",I4.4)'))
      endfor
  endelse

  return, result
end;function cmsystime
 

PRO polarimetric_processing_scan,scatvars,configvars,calvars,geometric_peak,range_peak_signal,scan_index,$
sweep_count,transition_flag,elevation,line_elevation,line_height,ngates,nlines,$
range_gate_spacing,range,height,gate_offset,gate_plot_max,spec,ainv,finv,l_matrix_vec,$
c_matrix_vec,l_chi_gamma,positioner_state,rho_hv_vec,phase_hv_deg_vec,total_power,range_centroid_signal,az_proc_index,sweep_count_override

;waitkey,dum
window,0,retain=2,xsize=550,ysize=300,xpos=140,ypos=125
window,1,retain=2,xsize=550,ysize=300,xpos=700,ypos=125
window,2,retain=2,xsize=550,ysize=300,xpos=140,ypos=450
window,3,retain=2,xsize=550,ysize=300,xpos=700,ypos=450

;
;finds peak of power-averaged range profiles (sweeps) 
;computes calibrated 2x2 scattering matrix, Smat, for all range gates
;within peak region
;computes 4x4 covariance matrix of Smat summed over peak region
;computes various derived products from covariance matrix
;scan_index=-5 when doing corner scan
;scan_index=-4 when staring  
;scan_index =-3 when looking at source horn
;scan_index=-2 when looking at sky
;scan_index =-1 when looking at corner reflector
;scan_index=1 when scanning ground



  nlines = max(sweep_count)+1   ; number of surface scan lines
  if nlines eq 0 then nlines=1
  start_line=0
  if configvars.i_el_override eq 1 then begin
     sweep_count_max=max(sweep_count)
     elvals=intarr(sweep_count_max+1)
    for i=0,sweep_count_max do begin
       index = where(sweep_count EQ i and transition_flag eq 0 and scan_index eq 1,count)
       elvals(i)=round(median(elevation(index)))
     endfor
    index0=where(elvals ge configvars.elmin_proc,count0)
    if count0 gt 0 then begin
       start_line=index0(0)
    endif else begin
        print,'elmin_proc is greater than all elevation angles'
        stop
     endelse
     index1=where(elvals gt configvars.elmax_proc,count1)
     if count1 gt 0 then stop_line=index1(0)-1 else begin
        print,'elmax_proc is greater than all elevation angles; processing up to max elevation'
        stop_line=sweep_count_max
     endelse
     nlines=stop_line-start_line+1

  endif
  if configvars.i_corner_process eq 1 then nlines=1 ; only look at corner reflector data if i_corner_process=1
  c_matrix_vec = complexarr(4,4,nlines) ;covariance matrix
  L_matrix_vec = fltarr(4,4,nlines)     ;Mueller matrix
  total_power=fltarr(nlines,4)          ;stores VV/HV/VH/HH power from ground target
  rho_hv_vec=fltarr(nlines)
  phase_hv_deg_vec=fltarr(nlines)
  rdepol_dB_vec=fltarr(nlines)
  line_index = lonarr(nlines)   ;sweep number at the beginning of lines having data (even numbered lines have data)
  block_count = lonarr(nlines)  ; number of VV/HV/VH/HH/cal/noise data blocks per sweep 
  line_elevation = fltarr(nlines) ;antenna elevation angle
  line_height=fltarr(nlines)
  range_peak_signal=fltarr(nlines)
  range_centroid_signal=fltarr(nlines)
  total_power=fltarr(nlines,4)  ;stores VV/HV/VH/HH power from ground target
  peak_power=fltarr(nlines,4)   ;stores VV/HV/VH/HH peak power from ground target

  if nlines gt 1 or configvars.i_corner_process eq 0 then begin 
;
;generate line_index, block_count (number of sweeps per azimuth cut),
;and line_elevation vectors (skips every other line to avoid short
;elevation segments between azimuth cuts)
;
 
     FOR iline0=0,nlines-1 DO BEGIN
        iline=iline0+start_line
        index = where(sweep_count EQ iline and transition_flag eq 0 and scan_index eq 1,count);added scan_index test 10/24/19
        if configvars.i_az_override eq 1 then begin
           index0=where(sweep_count_override eq iline,count)
           index=az_proc_index(index0);only process data over specified azimuth angle range
        endif
    if count eq 0 then begin
           count=1
           index(0)=n_elements(elevation)-1
        endif
        line_index(iline0) = index(0)
        block_count(iline0) =count
        line_elevation(iline0) = elevation(index(0))
        line_height(iline0)=height(index(0)) 
     ENDFOR 

  endif else begin
     index = where(scan_index eq -1,count)
     if count eq 0 then begin
        count=1
        index(0)=0
    endif
     line_index(0) = index(0)
     block_count(0) = 1         ;count
     line_elevation(0) = elevation(index(0))
  endelse
;
;loop through all lines to compute data products
;
smoothfac=configvars.smoothfac
if smoothfac lt 1 then smoothfac=1
if configvars.i_corner_process eq 1 and smoothfac ne 1 then begin
   smoothfac=1
   print,'setting smoothing factor to 1 when processing corner reflector data'
   waitkey,dum
endif

FOR iline=0,nlines-1 DO BEGIN
   n_blocks_per_line=block_count(iline)
   line_el=line_elevation(iline)
   process_polarimetric_data,scatvars,configvars,calvars,range,spec,n_blocks_per_line,ngates,smoothfac,c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_dB_vec,total_power,range_peak_signal,range_centroid_signal,iline,line_index,gate_offset,gate_plot_max,gate_peak,gate_peak_cal,ainv,finv,scan_index,line_el,peak_power
endfor;iline

wset,3
erase
;goto,skip_plot_1
wset,0
if nlines gt 1 then begin
   plot,line_elevation,rho_hv_vec,xtitle='elevation angle (deg)',ytitle='correlation coeff. mag.',title='VV/HH cor. coefficient vs. elevation angle',charsize=1.4  
wset,1
 plot,line_elevation,phase_hv_deg_vec,xtitle='elevation angle (deg)',ytitle='correlation coeff. phase (deg)',title='phase VV/HH cor. coe. vs elevation',charsize=1.4,yr=[-180,180]
wset,2
 plot,line_elevation,rdepol_dB_vec,xtitle='elevation angle (deg)',ytitle='depolarization ratio (dB)',title='depolarization vs elevation angle',charsize=1.4
waitkey,dum
endif

skip_plot_1:


wset,2
erase
;stop
END                             ;polarimetric_processing_scan procedure
   
PRO polarimetric_processing_stare,scatvars,configvars,calvars,ngates,range_gate_spacing,range,gate_offset,gate_plot_max,spec,ainv,finv,l_matrix_vec,$
c_matrix_vec,rho_hv_vec,phase_hv_deg_vec,total_power,n_blocks,elapsed_time,n_groups,range_peak_signal,range_centroid_signal,$
group_index,n_blocks_per_group,scan_index

;note: scan_index =-4 when staring


;waitkey,dum
scale=1.
window,0,retain=2,xsize=450*scale,ysize=300*scale,xpos=20,ypos=100
window,1,retain=2,xsize=450*scale,ysize=300*scale,xpos=500,ypos=100
window,2,retain=2,xsize=450*scale,ysize=300*scale,xpos=20,ypos=450
window,3,retain=2,xsize=450*scale,ysize=300*scale,xpos=500,ypos=450


;data is partitioned into groups of fixed number of data blocks
;determined by the group averaging time from the IDL configuration file.  This
;allows a large stare-mode file to be processed in sections of
;reasonable length in time

;
;two modes of processing are implemented:
;computing covariance matrix gate by gate  
;computing covariance matrix averaged over all gates


;gate-by-gate mode
;
;computes calibrated 2x2 scattering matrix, Smat, for all range gates
;within peak region
;
;computes 4x4 covariance matrix of Smat for all range gates (c_matrix_vs_bin)
;computes various derived products from covariance matrix for all range gates
;
;average mode
;
;averages covariance matrix over all gates to compute single covariance matrix (c_matrix)
;
n_bins_max=n_elements(range)
n_groups=ceil(max(elapsed_time)/configvars.group_averaging_time);eventually, couple this to time and/or distance tranvelled
n_blocks_per_group=ceil(n_blocks/n_groups)
c_matrix_vec = complexarr(4,4,n_groups) ;covariance matrix for average of all gates
L_matrix_vec = fltarr(4,4,n_groups) ;Mueller matrix for average of all gates
total_power=fltarr(n_groups,4)      ;stores VV/HV/VH/HH power from ground target
peak_power=fltarr(n_groups,4);stores VV/HV/VH/HH peak power from ground target
rho_hv_vec=fltarr(n_groups)
phase_hv_deg_vec=fltarr(n_groups)
rdepol_dB_vec=fltarr(n_groups)
range_peak_signal=fltarr(n_groups)
range_centroid_signal=fltarr(n_groups)

;
;loop through all groups to compute data products
;

group_index=lindgen(n_groups)*n_blocks_per_group

smoothfac=configvars.smoothfac

if smoothfac lt 1 then smoothfac=1
if configvars.i_corner_process eq 1 and smoothfac ne 1 then begin
   smoothfac=1
   print,'setting smoothing factor to 1 when processing corner reflector data'
   waitkey,dum
endif
line_el=0
FOR igroup=0,n_groups-1 DO BEGIN
   process_polarimetric_data,scatvars,configvars,calvars,range,spec,n_blocks_per_group,ngates,smoothfac,c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_db_vec,total_power,range_peak_signal,range_centroid_signal,igroup,group_index,gate_offset,gate_plot_max,gate_peak,gate_peak_cal,ainv,finv,scan_index,line_el,peak_power 
 endfor                         ;                         ;igroup

skip_plot_1:
wset,0
if n_groups gt 1 then begin
   plot,rho_hv_vec,xtitle='group number',ytitle='correlation coeff. mag.',title='VV/HH cor. coefficient vs. group',charsize=1.4  
   wset,1
   plot,phase_hv_deg_vec,xtitle='group number',ytitle='correlation coeff. phase (deg)',title='phase VV/HH cor. coe. vs group',charsize=1.4,yr=[-180,180]
   wset,2
   plot,rdepol_dB_vec,xtitle='group number',ytitle='depolarization ratio (dB)',title='depolarization vs group',charsize=1.4
   waitkey,dum
endif
goto,skip_env_chamber_plot
wset,2
erase
wset,3

peak_power_vv=20*alog10(abs(spec(*,0,gate_peak)))
peak_power_vh=20*alog10(abs(spec(*,1,gate_peak)))
peak_power_hv=20*alog10(abs(spec(*,2,gate_peak)))
peak_power_hh=20*alog10(abs(spec(*,3,gate_peak)))
peak_power_cal=20*alog10(abs(spec(*,4,gate_peak_cal)))
gate_max_noise=2*ngates-1
peak_power_noise=20*alog10(abs(spec(*,5,gate_max_noise)))
plot,elapsed_time/60,peak_power_vv,xtitle='elapsed time (min)',ytitle='dB',title='peak power vs time',charsize=1.5,yr=[-35,-5],xs=1,ys=1
oplot,elapsed_time/60,peak_power_cal,color=160
oplot,elapsed_time/60,peak_power_cal+mean(peak_power_vv)-mean(peak_power_cal),color=160,line=2
oplot,elapsed_time/60,peak_power_vh,color=40
oplot,elapsed_time/60,peak_power_hv,color=60
oplot,elapsed_time/60,peak_power_hh,color=80
oplot,elapsed_time/60,peak_power_noise-60,color=120

range_print=strmid(strtrim(range(gate_peak),2),0,4)+' m'
xyouts,.6*median(elapsed_time)/60,-8,'range: ',alignment=1,charsize=1.5
xyouts,.6*median(elapsed_time)/60,-8,range_print,alignment=0,charsize=1.5

cal_range_print=strmid(strtrim(range(gate_peak_cal),2),0,4)+' m'
xyouts,.9*median(elapsed_time)/60,-20,'cal loop range: ',alignment=1,charsize=1.5,color=160
xyouts,.9*median(elapsed_time)/60,-20,cal_range_print,alignment=0,charsize=1.5,color=160

wset,1
delta_p_vv_cal=peak_power_vv-peak_power_cal
delta_p_vh_cal=peak_power_vh-peak_power_cal
delta_p_hv_cal=peak_power_hv-peak_power_cal
delta_p_hh_cal=peak_power_hh-peak_power_cal

plot,elapsed_time/60,delta_p_vv_cal-mean(delta_p_vv_cal),xtitle='elapsed time (min)',ytitle='dB',title='corr. pwr: VV blk; HH red; VH blu; HV grn',charsize=1.5,yr=[-2,2],xs=1
oplot,elapsed_time/60,delta_p_vh_cal-mean(delta_p_vh_cal),color=40
oplot,elapsed_time/60,delta_p_hv_cal-mean(delta_p_hv_cal),color=60
oplot,elapsed_time/60,delta_p_hh_cal-mean(delta_p_hh_cal),color=160

stop
skip_env_chamber_plot:
erase


END                             ;polarimetric_processing_stare procedure
PRO polarimetric_processing_stare_by_independent_samples,scatvars,configvars,calvars,ngates,range_gate_spacing,range,$
   gate_offset,gate_plot_max,spec,ainv,finv,l_matrix_vec,c_matrix_vec,rho_hv_vec,phase_hv_deg_vec,total_power,n_blocks,$
   elapsed_time,n_groups,range_peak_signal,range_centroid_signal,group_index,n_blocks_per_group,scan_index,$
   independent_sample_index,distance

  
;note: scan_index =-4 when staring


;waitkey,dum
scale=1.
window,0,retain=2,xsize=450*scale,ysize=300*scale,xpos=20,ypos=100
window,1,retain=2,xsize=450*scale,ysize=300*scale,xpos=500,ypos=100
window,2,retain=2,xsize=450*scale,ysize=300*scale,xpos=20,ypos=450
window,3,retain=2,xsize=450*scale,ysize=300*scale,xpos=500,ypos=450


;data is partitioned into groups of fixed number of data blocks
;determined by the group averaging time from the IDL configuration file.  This
;allows a large stare-mode file to be processed in sections of
;reasonable length in time

;
;two modes of processing are implemented:
;computing covariance matrix gate by gate  
;computing covariance matrix averaged over all gates
 

;gate-by-gate mode
;
;computes calibrated 2x2 scattering matrix, Smat, for all range gates
;within peak region
;
;computes 4x4 covariance matrix of Smat for all range gates (c_matrix_vs_bin)
;computes various derived products from covariance matrix for all range gates
;
;average mode
;
;averages covariance matrix over all gates to compute single covariance matrix (c_matrix)
;

n_blocks_ind=n_elements(independent_sample_index)
n_bins_max=n_elements(range)
n_groups=floor(n_blocks_ind/configvars.n_blocks_ind_per_group)
n_blocks_per_group=configvars.n_blocks_ind_per_group
c_matrix_vec = complexarr(4,4,n_groups) ;covariance matrix for average of all gates
L_matrix_vec = fltarr(4,4,n_groups) ;Mueller matrix for average of all gates
total_power=fltarr(n_groups,4)      ;stores VV/HV/VH/HH power from ground target
peak_power=fltarr(n_groups,4);stores VV/HV/VH/HH peak power from ground target
rho_hv_vec=fltarr(n_groups)
phase_hv_deg_vec=fltarr(n_groups)
rdepol_dB_vec=fltarr(n_groups)
range_peak_signal=fltarr(n_groups)
range_centroid_signal=fltarr(n_groups)
hh_image0=fltarr(n_blocks_ind,ngates)
vh_image0=fltarr(n_blocks_ind,ngates)
hv_image0=fltarr(n_blocks_ind,ngates)
vv_image0=fltarr(n_blocks_ind,ngates)
;
;loop through all groups to compute data products
;

group_index=lindgen(n_groups)*n_blocks_per_group

smoothfac=configvars.smoothfac

if smoothfac lt 1 then smoothfac=1
if configvars.i_corner_process eq 1 and smoothfac ne 1 then begin
   smoothfac=1
   print,'setting smoothing factor to 1 when processing corner reflector data'
   waitkey,dum
endif

;create matrix of spectra that are independent

spec_independent=spec(independent_sample_index,*,*)
line_el=0
FOR igroup=0,n_groups-1 DO BEGIN
   
   wset,2
   ;highlight along track transect position of group being processed
   plot,elapsed_time/60,distance,charsize=1.5,xtitle='time (min)',ytitle='distance (m)'
   samplevec=independent_sample_index(group_index(igroup):group_index(igroup)+n_blocks_per_group)
   oplot,elapsed_time(independent_sample_index)/60,distance(independent_sample_index),psym=1,color=60 
   oplot,elapsed_time(samplevec)/60,distance(samplevec),psym=1,color=160
  

   process_polarimetric_data,scatvars,configvars,calvars,range,spec_independent,n_blocks_per_group,ngates,smoothfac,c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_db_vec,total_power,range_peak_signal,range_centroid_signal,igroup,group_index,gate_offset,gate_plot_max,gate_peak,gate_peak_cal,ainv,finv,scan_index,line_el,peak_power 
   HH_image0(group_index(igroup):group_index(igroup)+n_blocks_per_group,*)=20*alog10(abs(spec_independent(group_index(igroup):group_index(igroup)+n_blocks_per_group,3,0:ngates-1)))
   VV_image0(group_index(igroup):group_index(igroup)+n_blocks_per_group,*)=20*alog10(abs(spec_independent(group_index(igroup):group_index(igroup)+n_blocks_per_group,0,0:ngates-1)))
   HV_image0(group_index(igroup):group_index(igroup)+n_blocks_per_group,*)=20*alog10(abs(spec_independent(group_index(igroup):group_index(igroup)+n_blocks_per_group,1,0:ngates-1)))
  VH_image0(group_index(igroup):group_index(igroup)+n_blocks_per_group,*)=20*alog10(abs(spec_independent(group_index(igroup):group_index(igroup)+n_blocks_per_group,2,0:ngates-1)))
endfor                          ;                         ;igroup

;trim off zeros
HH_image=HH_image0(0:max(group_index)+n_blocks_per_group-1,*)
VV_image=VV_image0(0:max(group_index)+n_blocks_per_group-1,*)
HV_image=HV_image0(0:max(group_index)+n_blocks_per_group-1,*)
VH_image=VH_image0(0:max(group_index)+n_blocks_per_group-1,*)
LDR_image=(VH_image+HV_image)/2-(HH_image+VV_image)/2
HHVV_ave_image=(HH_image+VV_image)/2

;generate 2D plots of scattered power
window,1,re=2,xsize=600,ysize=450,xpos=100,ypos=100
min_range=1.0;m
max_range=3.0;m
index_range_min=where(range gt min_range)
index_range_max=where(range gt max_range)
;set up plotting parameters
x_range=[min(distance(independent_sample_index)),max(distance(independent_sample_index))]
y_range=[range(index_range_max(0)),range(index_range_min(0))]

xlabel='approximate distance along transect (m)'
ylabel='radar range from antenna (m)'
smoothfac=1
replot:
hh_image_smooth=congrid(smooth(hh_image(*,index_range_min(0):index_range_max(0)),[smoothfac,1],/edge_truncate),200,200)
hv_image_smooth=congrid(smooth(hv_image(*,index_range_min(0):index_range_max(0)),[smoothfac,1],/edge_truncate),200,200)
vh_image_smooth=congrid(smooth(vh_image(*,index_range_min(0):index_range_max(0)),[smoothfac,1],/edge_truncate),200,200)
vv_image_smooth=congrid(smooth(vv_image(*,index_range_min(0):index_range_max(0)),[smoothfac,1],/edge_truncate),200,200)
hhvv_ave_image_smooth=congrid(smooth(hhvv_ave_image(*,index_range_min(0):index_range_max(0)),[smoothfac,1],/edge_truncate),200,200)
LDR_image_smooth=congrid(smooth(LDR_image(*,index_range_min(0):index_range_max(0)),[smoothfac,1],/edge_truncate),200,200)

data_max=max(hh_image)
data_min=data_max-50.0;dB
ct_string='relative power (dB)'

vplot=reverse(hh_image_smooth,2)
title_str='HH power along transect'
erase
image_plot,vplot,data_min,data_max,x_range,y_range,ct_string,xlabel,ylabel,title_str
waitkey,dum

vplot=reverse(vv_image_smooth,2)
title_str='VV power along transect'
erase
image_plot,vplot,data_min,data_max,x_range,y_range,ct_string,xlabel,ylabel,title_str
waitkey,dum

;vplot=reverse(hhvv_ave_image_smooth,2)
;title_str='HH/VV average power along transect'
;erase
;image_plot,vplot,data_min,data_max,x_range,y_range,ct_string,xlabel,ylabel,title_str
;waitkey,dum

ct_string='LDR (dB)'
data_max=5
data_min=-40
vplot=reverse(LDR_image_smooth,2)
title_str='Linear Depolarization Ratio along transect'
erase
image_plot,vplot,data_min,data_max,x_range,y_range,ct_string,xlabel,ylabel,title_str
waitkey,dum

if configvars.i_batch eq 0 then read,'enter smoothing factor for 2D image, 99 to exit: ',smoothfac else smoothfac=99

if smoothfac ne 99 then goto,replot
window,1,retain=2,xsize=450*scale,ysize=300*scale,xpos=500,ypos=100
wset,0
wshow,0
wshow,2
wshow,3
if n_groups gt 1 then begin
   plot,rho_hv_vec,xtitle='group number',ytitle='correlation coeff. mag.',title='VV/HH cor. coefficient vs. group',charsize=1.4  
   wset,1
   plot,phase_hv_deg_vec,xtitle='group number',ytitle='correlation coeff. phase (deg)',title='phase VV/HH cor. coe. vs group',charsize=1.4,yr=[-180,180]
   wset,2
   plot,rdepol_dB_vec,xtitle='group number',ytitle='depolarization ratio (dB)',title='depolarization vs group',charsize=1.4
   waitkey,dum
endif
for i=0,3 do begin
   wset,i
   erase
endfor


END                             ;polarimetric_processing_stare_by_independent_samples procedure
PRO nrcs_compute_scan,scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,range_centroid_signal,corner_range,pos_height,$
line_elevation,nlines,reference_calibration_loop_power,current_calibration_loop_power,L_matrix,processed_data_filename,rho_hv_vec,$
phase_hv_deg_vec,total_power,line_height,total_corner_power_vv,total_corner_power_hh,gps_latitude,gps_longitude,calfile

for i=0,3 do begin
   wdelete,i
endfor

total_corner_power_vv_file=10^(calvars.corner_reflector_vv_power_dBm/10.)
total_corner_power_hh_file=10^(calvars.corner_reflector_hh_power_dBm/10.)
reference_calibration_loop_power_file=10^(calvars.cal_peak_dBm/10.)
corner_range_file=calvars.corner_reflector_range_m

if configvars.i_corner_cal_override eq 0 then begin
  total_corner_power_vv_in_use=total_corner_power_vv_file
  total_corner_power_hh_in_use=total_corner_power_hh_file
  reference_calibration_loop_power_in_use=reference_calibration_loop_power_file
  corner_range_in_use=corner_range_file
endif
if configvars.i_corner_cal_override eq 1 then begin
  total_corner_power_vv_in_use=total_corner_power_vv
  total_corner_power_hh_in_use=total_corner_power_hh
  reference_calibration_loop_power_in_use=reference_calibration_loop_power
  corner_range_in_use=corner_range
endif


if nlines gt 1 then begin
   window,0,retain=2,xsize=550,ysize=300,xpos=140,ypos=125
   window,1,retain=2,xsize=550,ysize=300,xpos=700,ypos=125
   window,2,retain=2,xsize=550,ysize=300,xpos=140,ypos=450
endif

;window,3,retain=2,xsize=550,ysize=300,xpos=700,ypos=450
;window,4,retain=2,xsize=550,ysize=215,xpos=700,ypos=0

;computes normalized radar cross section of target and corrects for near-field
;gain loss

   pi = !pi
   ln_2 = alog10(2)/alog10(exp(1))
   theta_i = line_elevation*!dtor
   nrcs = fltarr(4,nlines)
      nrcs_alt = fltarr(4,nlines)
   antenna_beamwidth_rad = private_config.antenna_beam_width*!dtor

   corr_cal = current_calibration_loop_power/reference_calibration_loop_power_in_use
   if configvars.i_cal_loop_override eq 0 then corr_cal=1.0

   print,'******************************'
   print,'reference calibration power from data file (dB):',10*alog10(reference_calibration_loop_power_file)
   print,'reference calibration power from cal file(dB):',10*alog10(reference_calibration_loop_power)  
   print,'reference calibration power used for NRCS calculation (dB):',10*alog10(reference_calibration_loop_power_in_use)  
   print,'current calibration power (dB): ',10*alog10(current_calibration_loop_power)
   if configvars.i_cal_loop_override eq 0 then begin
      print,'drift correction factor (linear): ',corr_cal
   endif else begin
      print,'i_cal_loop_override=1, no drift correction applied'
   endelse
   print,'*****************************'
   scale_factor = fltarr(nlines)
   scale_factor_alt = fltarr(nlines)
   corner_sigma=private_config.corner_reflector_sigma

   FOR iline=0,nlines-1 DO BEGIN 
;
;using range centroid instead of h/cos(theta_i) as a more accurate
;estimate of range
;      scale_factor_alt(iline) = 8*ln_2*line_height(iline)^2*corner_sigma/(pi*corner_range^4*antenna_beamwidth_rad^2*cos(theta_i(iline)))/corr_cal; 
      scale_factor(iline) = 8*ln_2*range_centroid_signal(iline)^2*corner_sigma*cos(theta_i(iline))/(pi*corner_range_in_use^4*antenna_beamwidth_rad^2)/corr_cal; 

      if configvars.i_corner_process eq 1 then scale_factor(iline)=range_peak_signal(iline)^4*corner_sigma/corner_range_file^4

      
      FOR j=0,3 DO BEGIN 
         
         if j eq 0 then cr_power=total_corner_power_vv_in_use

         if j eq 1 or j eq 2 then cr_power=(total_corner_power_vv_in_use+total_corner_power_hh_in_use)/2

        if j eq 3 then cr_power=total_corner_power_hh_in_use
 
   pratio =total_power(iline,j)/cr_power ;float(c_matrix(j,j,iline))
   nrcs(j,iline)= scale_factor(iline)*pratio
;   nrcs_alt(j,iline)=scale_factor_alt(iline)*pratio
   
      ENDFOR
  ENDFOR 
;average cross-pol data (we know that NRCS of cross-pol terms should
;be equal)


nrcs_cross_pol=(nrcs(1,*)+nrcs(2,*))/2.0
;nrcs(1,*)=nrcs_cross_pol
;nrcs(2,*)=nrcs_cross_pol 
nrcs_db=10*alog10(nrcs)
;nrcs_alt_db=10*alog10(nrcs_alt)
if nlines gt 1 then begin
   wset,0
plot,line_elevation(*),nrcs_dB(0,*),yrange=[-60,10],xtitle='elevation angle, deg',ytitle='NRCS, dB m2/m2',charsize=1.6,title='NRCS VV blk; HH red; VH blue; HV grn;',xs=1
   oplot,line_elevation(*),nrcs_dB(1,*),color=80,line=2
   oplot,line_elevation(*),nrcs_dB(2,*),color=40,line=2
   oplot,line_elevation(*),nrcs_dB(3,*),color=160
;   oplot,line_elevation(*),nrcs_alt_dB(0,*),color=25,line=1,thick=2
   wset,1
   plot,line_elevation,rho_hv_vec,xtitle='elevation angle (deg)',ytitle='correlation coeff. mag.',title='VV/HH cor. coefficient vs. elevation angle',charsize=1.4  
wset,2
 plot,line_elevation,phase_hv_deg_vec,xtitle='elevation angle (deg)',ytitle='correlation coeff. phase (deg)',title='phase VV/HH cor. coe. vs elevation',charsize=1.4,yr=[-180,180],ys=1
endif 
;store processed data
extension='.scan'
if configvars.i_corner_process eq 1 then extension='.scan.corner'
   print,'summary data stored in: ',processed_data_filename+extension
   openw,1,processed_data_filename+extension
   printf,1,'number of elevation angles: '
   printf,1,nlines
   printf,1
   printf,1,'elevation angles (deg)'
   printf,1,line_elevation
   printf,1
   printf,1,'range peak signal (m): '
   printf,1,range_peak_signal
   printf,1
   printf,1,'corner reflector sigma (square meters): '
   printf,1,corner_sigma
   printf,1
   printf,1,'6 dB two-way antenna beamwidth (deg): '
   printf,1, private_config.antenna_beam_width
   printf,1
   printf,1,'current cal loop power (dB): '
   printf,1,10*alog10(current_calibration_loop_power)
   printf,1
   printf,1,'reference cal loop power (dB): '
   printf,1,10*alog10(reference_calibration_loop_power)
   printf,1
   printf,1,'corner range (m): '
   printf,1,corner_range
   printf,1
   printf,1,'latitude: '
   printf,1,mean(gps_latitude)
   printf,1
   printf,1,'longitude: '
   printf,1,mean(gps_longitude)
   printf,1
   printf,1,'Mueller matrix for various elevation angles: '
   for iline=0,nlines-1 do begin
      for i=0,3 do begin      
         L_matrix_norm=L_matrix(*,*,iline)/max(L_matrix(*,*,iline))
         printf,1,L_matrix_norm(0,i),' ',L_matrix_norm(1,i),' ',L_matrix_norm(2,i),' ',L_matrix_norm(3,i),' '
      endfor
      printf,1
   endfor
   
   printf,1,'Covariance matrix for various elevation angles: '
   for iline=0,nlines-1 do begin
      for i=0,3 do begin      
        c_matrix_norm=c_matrix(*,*,iline)/max(abs(c_matrix(*,*,iline)))
        real0=real_part(c_matrix_norm(0,i))
        real1=real_part(c_matrix_norm(1,i))
        real2=real_part(c_matrix_norm(2,i))
        real3=real_part(c_matrix_norm(3,i))
        imag0=imaginary(c_matrix_norm(0,i))
        imag1=imaginary(c_matrix_norm(1,i))
        imag2=imaginary(c_matrix_norm(2,i))
        imag3=imaginary(c_matrix_norm(3,i))
        printf,1,'('+strtrim(string(real0),2)+', '+strtrim(string(imag0),2)+') ('+strtrim(string(real1),2)+',  '+strtrim(string(imag1),2)+') ('+strtrim(string(real2),2)+', '+strtrim(string(imag2),2)+') ('+strtrim(string(real3),2)+', '+strtrim(string(imag3),2)+')'
        

      endfor
      printf,1
   endfor
   
   printf,1,'magnitude of HH/VV correlation coefficient: '
   printf,1,rho_hv_vec
   printf,1
   printf,1,'phase of HH/VV correlation coefficient (deg): '
   printf,1,phase_hv_deg_vec
   printf,1
   if configvars.i_corner_process eq 0 then printf,1,'Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH; rows: elevations angles : '
   if configvars.i_corner_process eq 1 then printf,1,'Radar Cross Section (dBm2): VV, HV, VH, HH: '
   printf,1,nrcs_db                     ; 
   printf,1                             ; 
   printf,1,'calibration filename'
   printf,1,calfile
   close,1
   print,'done'
;waitkey,dum

END                             ;nrcs_compute_scan procedure

PRO  nrcs_compute_stare,scatvars,calvars,configvars,private_config,c_matrix,range_peak_signal,range_centroid_signal,$
corner_range,pos_height,n_groups,reference_calibration_loop_power,current_calibration_loop_power,L_matrix,processed_data_filename,rho_hv_vec,$
phase_hv_deg_vec,total_power,total_corner_power_vv,total_corner_power_hh,calfile,time_sec,gps_latitude,gps_longitude,group_index,along_track_tilt,cross_track_tilt,n_blocks_per_group,independent_sample_index
for i=0,3 do begin
wdelete,i
endfor


if n_groups gt 1 then begin
   window,0,retain=2,xsize=550,ysize=300,xpos=140,ypos=125
   window,1,retain=2,xsize=550,ysize=300,xpos=700,ypos=125
   window,2,retain=2,xsize=550,ysize=300,xpos=140,ypos=450
endif

;
;check if stored calibration file has same corner power as current file

total_corner_power_vv_file=10^(calvars.corner_reflector_vv_power_dBm/10.)
total_corner_power_hh_file=10^(calvars.corner_reflector_hh_power_dBm/10.)
corner_range_file=calvars.corner_reflector_range_m
reference_calibration_loop_power_file=10^(calvars.cal_peak_dBm/10.)

goto,jump_message
if abs(calvars.corner_reflector_vv_power_dBm-10*alog10(total_corner_power_vv)) gt 0.1 or abs(calvars.corner_reflector_hh_power_dBm-10*alog10(total_corner_power_hh)) gt 0.1 then begin

print 
print,'corner powers differ between current calibration file and values in data file: '
print,'data file values: '
print,'corner reflector power vv: ',total_corner_power_vv_file 
print,'corner reflector power hh: ',total_corner_power_hh_file 
print,'corner reflector range: ',corner_range_file
print
print,'calibration file values: '
print,'corner reflector power vv: ',total_corner_power_vv 
print,'corner reflector power hh: ',total_corner_power_hh 
print,'corner reflector range: ',corner_range
print
print,'this is just a heads up.'

waitkey,dum
endif
jump_message:
;computes normalized radar cross section of target 

   pi = !pi
   ln_2 = alog10(2)/alog10(exp(1))

   nrcs = fltarr(4,n_groups)
   antenna_beamwidth_rad = private_config.antenna_beam_width*!dtor

   corr_cal = current_calibration_loop_power/reference_calibration_loop_power_file
   print,'******************************'
   print,'reference calibration power from data file (dB):',10*alog10(reference_calibration_loop_power_file)
   print,'reference calibration power from cal file(dB):',10*alog10(reference_calibration_loop_power)  
   print,'current calibration power (dB): ',10*alog10(current_calibration_loop_power)
   print,'drift correction factor (linear): ',corr_cal
   print,'*****************************'

   corner_sigma=private_config.corner_reflector_sigma
   scale_factor=fltarr(n_groups)

                                ;find mean tilt of instrument; psi is
                                ;angle from z-axis to point on surface
   psi=atan(sqrt(tan(along_track_tilt*!dtor)^2+tan(cross_track_tilt*!dtor)^2))


   FOR igroup=0,n_groups-1 DO BEGIN 
      vec=group_index(igroup)+indgen(n_blocks_per_group)   
      mean_tilt=mean(psi(vec))
      scale_factor(igroup) = 8*ln_2*range_centroid_signal(igroup)^2*corner_sigma*cos(mean_tilt)/(pi*corner_range_file^4*antenna_beamwidth_rad^2)/corr_cal ; 
      if configvars.i_corner_process eq 1 then scale_factor(igroup)=range_peak_signal(igroup)^4*corner_sigma/corner_range_file^4      
      FOR j=0,3 DO BEGIN 

         if j eq 0 then cr_power=total_corner_power_vv_file

         if j eq 1 or j eq 2 then cr_power=(total_corner_power_vv_file+total_corner_power_hh_file)/2

        if j eq 3 then cr_power=total_corner_power_hh_file
 
        pratio =total_power(igroup,j)/cr_power
      
         nrcs(j,igroup)= scale_factor(igroup)*pratio
      ENDFOR
  ENDFOR 




nrcs_db=10*alog10(nrcs)
if n_groups gt 1 then begin
   wset,0
plot,nrcs_dB(0,*),yrange=[-60,10],xtitle='group number',ytitle='NRCS, dB m2/m2',charsize=1.6,title='NRCS VV blk; HH red; VH blue; HV grn;',xs=1
   oplot,nrcs_dB(1,*),color=80,line=2
   oplot,nrcs_dB(2,*),color=40,line=2
   oplot,nrcs_dB(3,*),color=160
wset,1
   plot,rho_hv_vec,xtitle='group number',ytitle='correlation coeff. mag.',title='VV/HH cor. coefficient vs. group',charsize=1.4  
wset,2
 plot,phase_hv_deg_vec,xtitle='group number',ytitle='correlation coeff. phase (deg)',title='phase VV/HH cor. coe. vs group',charsize=1.4,yr=[-180,180],ys=1
endif else begin
print
print,'NRCS (dB): '
print,nrcs_dB
print,'scale factor: ',scale_factor

endelse

;store processed data
extension='.stare'
if configvars.i_corner_process eq 1 then extension='stare.corner'
   print,'summary data stored in: ',processed_data_filename+extension
   openw,1,processed_data_filename+extension
   printf,1,'number of data groups: '
   printf,1,n_groups
   printf,1
   printf,1,'n_blocks_per_group'
   printf,1,n_blocks_per_group
   printf,1
   printf,1,'range peak signal: '
   printf,1,range_peak_signal
   printf,1
   printf,1,'range centroid signal: '
   printf,1,range_centroid_signal
   printf,1
   printf,1,'corner reflector sigma (square meters): '
   printf,1,corner_sigma
   printf,1
   printf,1,'6 dB two-way antenna beamwidth (deg): '
   printf,1, private_config.antenna_beam_width
   printf,1
   printf,1,'current cal loop power (dB): '
   printf,1,10*alog10(current_calibration_loop_power)
   printf,1
   printf,1,'reference cal loop power (dB): '
   printf,1,10*alog10(reference_calibration_loop_power)
   printf,1
   printf,1,'corner range (m): '
   printf,1,corner_range
   printf,1
   for igroup=0,n_groups-1 do begin
      vec0=group_index(igroup)+indgen(n_blocks_per_group) 
      vec=independent_sample_index(vec0)
         printf,1,'group number: '
         printf,1,igroup
         printf,1
         printf,1,'group index: '
         printf,1,group_index(igroup)
         printf,1
         printf,1,'start time of group: '
         printf,1,time_sec(group_index(igroup))
         printf,1
         printf,1,'end time of group: '
         printf,1,time_sec(group_index(igroup)+n_blocks_per_group-1)
         printf,1
         printf,1,'mean latitude: '
         printf,1,mean(gps_latitude(vec))
         printf,1
         printf,1,'mean longitude: '
         printf,1,mean(gps_longitude(vec))
         printf,1
         printf,1,'mean along track tilt: '
         printf,1,mean(along_track_tilt(vec))
         printf,1
         printf,1,'standard deviation along track tilt: '
         printf,1,stddev(along_track_tilt(vec))
         printf,1
         printf,1,'mean cross track tilt: '
         printf,1,mean(cross_track_tilt(vec))
         printf,1
         printf,1,'standard deviation along track tilt: '
         printf,1,stddev(cross_track_tilt(vec))
         printf,1
         printf,1,'mean tilt relative to vertical axis: '
         printf,1,mean(psi(vec)/!dtor)
         printf,1
         printf,1,'standard deviation of tilt relative to vertical axis: '
         printf,1,stddev(psi(vec)/!dtor)
         printf,1
         printf,1,'Mueller matrix '
      for i=0,3 do begin  
         L_matrix_norm=L_matrix(*,*,igroup)/max(L_matrix(*,*,igroup))
         printf,1,L_matrix_norm(0,i),' ',L_matrix_norm(1,i),' ',L_matrix_norm(2,i),' ',L_matrix_norm(3,i),' '
      endfor
      printf,1
      printf,1,'normalized covariance matrix'
      for i=0,3 do begin      
         c_matrix_norm=c_matrix(*,*,igroup)/max(abs(c_matrix(*,*,igroup)))
        real0=real_part(c_matrix_norm(0,i))
        real1=real_part(c_matrix_norm(1,i))
        real2=real_part(c_matrix_norm(2,i))
        real3=real_part(c_matrix_norm(3,i))
        imag0=imaginary(c_matrix_norm(0,i))
        imag1=imaginary(c_matrix_norm(1,i))
        imag2=imaginary(c_matrix_norm(2,i))
        imag3=imaginary(c_matrix_norm(3,i))
        printf,1,'('+strtrim(string(real0),2)+', '+strtrim(string(imag0),2)+') ('+strtrim(string(real1),2)+',  '+strtrim(string(imag1),2)+') ('+strtrim(string(real2),2)+', '+strtrim(string(imag2),2)+') ('+strtrim(string(real3),2)+', '+strtrim(string(imag3),2)+')'
        
     endfor
      printf,1
      printf,1,'magnitude of HH/VV correlation coefficient: '
      printf,1,rho_hv_vec(igroup)
      printf,1
      printf,1,'phase of HH/VV correlation coefficient (deg): '
      printf,1,phase_hv_deg_vec(igroup)
      printf,1
      if configvars.i_corner_process eq 0 then printf,1,'Normalized Radar Cross Section (dBm2/m2) columns: VV, HV, VH, HH; rows: elevations angles : '
      if configvars.i_corner_process eq 1 then printf,1,'Radar Cross Section (dBm2): VV, HV, VH, HH: '
      printf,1,nrcs_db(*,igroup) ; 
      printf,1                   ;
   endfor                        ;
   printf,1,'calibration filename'
   printf,1,calfile
   close,1
   print,'done'

;waitkey,dum
END                             ;nrcs_compute_stare procedure

PRO read_configuration,configvars,iscat
print,'reading SCAT_IDL_config.txt'
configvars={instrument_name:strarr(1),$
i_calibrate:0B,$
calfile:strarr(1),$
i_corner_cal_override:0B,$
i_cal_loop_override:0B,$
i_batch:0B,$
i_raw_plot:0B,$
max_display_range: 0.0,$
i_proc_ind:0B,$
n_blocks_ind_per_group:0S,$
i_az_override: 0B,$
azmin_proc: 0.0,$
azmax_proc: 0.0,$
i_el_override: 0B,$
elmin_proc: 0.0,$
elmax_proc: 0.0,$
i_corner_process:0B,$
i_pol_scat:0B,$
i_pol_signature:0B,$
i_pol_data_in_footprint:0B,$
smoothfac:0B,$
i_temp_time_plot:0B,$
i_sub_band:0B,$
proc_thresh_left_dB:0.0,$
proc_thresh_right_dB:0.0,$
group_averaging_time:0.0,$
sub_bandwidth:0.0,$
sub_center_frequency:0.0,$   
raw_data_path: strarr(1),$
processed_data_path: strarr(1)$
}

;  read,'enter 1 for Ku-Scat, 2 for Ka-Scat: ',iscat
  
  if iscat eq 1 then filename='Ku-SCAT_IDL_config.txt'
  if iscat eq 2 then filename='Ka-SCAT_IDL_config.txt'

  stringvar=''
  bytevar=0B
  intvar=0S
  longvar=0L
  floatvar=0.0
  title=''
  close,1
  openr,1,filename
  readf,1,title

  test_string='instrument_name'
  get_variable,filename,stringvar,test_string
  configvars.instrument_name=stringvar

if iscat eq 1 and configvars.instrument_name ne 'Ku-Scat' then begin
   print
   print,'WARNING: instrument name in Ku-SCAT_IDL_config.txt must be Ku-Scat'
   print
   stop
end
if iscat eq 2 and configvars.instrument_name ne 'Ka-Scat' then begin
   print
   print,'WARNING: instrument name in Ka-SCAT_IDL_config.txt must be Ka-Scat'
   print
   stop
end

  test_string='i_sub_band'
  get_variable,filename,bytevar,test_string
  configvars.i_sub_band=bytevar

  test_string='sub_bandwidth'
  get_variable,filename,floatvar,test_string
  configvars.sub_bandwidth=floatvar

  test_string='sub_center_frequency'
  get_variable,filename,floatvar,test_string
  configvars.sub_center_frequency=floatvar

  test_string='i_calibrate'
  get_variable,filename,bytevar,test_string
  configvars.i_calibrate=bytevar

  test_string='calfile'
  get_variable,filename,stringvar,test_string
  configvars.calfile=stringvar

  test_string='i_corner_cal_override'
  get_variable,filename,bytevar,test_string
  configvars.i_corner_cal_override=bytevar

  test_string='i_cal_loop_override'
  get_variable,filename,bytevar,test_string
  configvars.i_cal_loop_override=bytevar

  test_string='i_batch'
  get_variable,filename,bytevar,test_string
  configvars.i_batch=bytevar

  test_string='i_raw_plot'
  get_variable,filename,bytevar,test_string
  configvars.i_raw_plot=bytevar

  test_string='max_display_range'
  get_variable,filename,floatvar,test_string
  configvars.max_display_range=floatvar

  test_string='i_proc_ind'
  get_variable,filename,bytevar,test_string
  configvars.i_proc_ind=bytevar

  test_string='n_blocks_ind_per_group'
  get_variable,filename,intvar,test_string
  configvars.n_blocks_ind_per_group=intvar

  test_string='i_az_override'
  get_variable,filename,bytevar,test_string
  configvars.i_az_override=bytevar

  test_string='azmin_proc'
  get_variable,filename,floatvar,test_string
  configvars.azmin_proc=floatvar

  test_string='azmax_proc'
  get_variable,filename,floatvar,test_string
  configvars.azmax_proc=floatvar

  test_string='i_el_override'
  get_variable,filename,bytevar,test_string
  configvars.i_el_override=bytevar

  test_string='elmin_proc'
  get_variable,filename,floatvar,test_string
  configvars.elmin_proc=floatvar

  test_string='elmax_proc'
  get_variable,filename,floatvar,test_string
  configvars.elmax_proc=floatvar

  test_string='i_corner_process'
  get_variable,filename,bytevar,test_string
  configvars.i_corner_process=bytevar

  test_string='i_pol_scat'
  get_variable,filename,bytevar,test_string
  configvars.i_pol_scat=bytevar

  test_string='i_pol_signature'
  get_variable,filename,bytevar,test_string
  configvars.i_pol_signature=bytevar
  
  test_string='i_pol_data_in_footprint'
  get_variable,filename,bytevar,test_string
  configvars.i_pol_data_in_footprint=bytevar

  test_string='smoothfac'
  get_variable,filename,bytevar,test_string
  configvars.smoothfac=bytevar
  
  test_string='i_temp_time_plot'
  get_variable,filename,bytevar,test_string
  configvars.i_temp_time_plot=bytevar

  test_string='proc_thresh_left_dB'
  get_variable,filename,floatvar,test_string
  configvars.proc_thresh_left_dB=floatvar

  test_string='proc_thresh_right_dB'
  get_variable,filename,floatvar,test_string
  configvars.proc_thresh_right_dB=floatvar

  test_string='group_averaging_time'
  get_variable,filename,floatvar,test_string
  configvars.group_averaging_time=floatvar

  test_string='raw_data_path'
  get_variable,filename,stringvar,test_string
  configvars.raw_data_path=stringvar

  test_string='processed_data_path'
  get_variable,filename,stringvar,test_string
  configvars.processed_data_path=stringvar

  close,1

  end;read_configuration   

PRO mueller_compute,c_matrix,c_matrix_vs_bin,L_matrix,L_matrix_vs_bin

  
;find L_matrix averaged over all ranges from individual s-matrix covariances

;note on notation: s_hhs=s_hh*
    s_vv2 = abs(c_matrix(0,0))
    s_vv_s_vhs =c_matrix(1,0)
    s_vv_s_hvs =c_matrix(2,0)
    s_vv_s_hhs = c_matrix(3,0)

    s_vh_s_vvs = c_matrix(0,1)
    s_vh2 =  abs(c_matrix(1,1))
    s_vh_s_hvs = c_matrix(2,1)
    s_vh_s_hhs = c_matrix(3,1)

    s_hv_s_vvs = c_matrix(0,2)
    s_hv_s_vhs =c_matrix(1,2)    
    s_hv2 = abs(c_matrix(2,2))
    s_hv_s_hhs =c_matrix(3,2)
    
    s_hh_s_vvs = c_matrix(0,3)
    s_hh_s_vhs =c_matrix(1,3)    
    s_hh_s_hvs = c_matrix(2,3)
    s_hh2 = abs(c_matrix(3,3))        

    L_matrix(0,0) = .5*(s_vv2+s_hh2+s_hv2+s_vh2)
    L_matrix(1,0) = .5*(s_vv2-s_hh2+s_hv2-s_vh2)
    L_matrix(2,0) = real_part(s_hv_s_hhs)+real_part(s_vv_s_vhs)
    L_matrix(3,0) = imaginary(s_hh_s_hvs)+imaginary(s_vh_s_vvs)
    
    L_matrix(0,1) = .5*(s_vv2-s_hh2-s_hv2+s_vh2)
    L_matrix(1,1) = .5*(s_vv2+s_hh2-s_hv2-s_vh2)
    L_matrix(2,1) = real_part(s_vv_s_vhs)-real_part(s_hh_s_hvs)
    L_matrix(3,1) = imaginary(s_vh_s_vvs)+imaginary(s_hv_s_hhs)
    
    L_matrix(0,2) = real_part(s_vh_s_hhs)+real_part(s_vv_s_hvs)
    L_matrix(1,2) = real_part(s_vv_s_hvs)-real_part(s_hh_s_vhs)
    L_matrix(2,2) = real_part(s_vv_s_hhs)+real_part(s_vh_s_hvs)
    L_matrix(3,2) = imaginary(s_hh_s_vvs)+imaginary(s_vh_s_hvs)
    
    L_matrix(0,3) = imaginary(s_vh_s_hhs)+imaginary(s_vv_s_hvs)
    L_matrix(1,3) = imaginary(s_vv_s_hvs)-imaginary(s_vh_s_hhs)
    L_matrix(2,3) = -(imaginary(s_hh_s_vvs)+imaginary(s_hv_s_vhs));negative sign added 8/1/19 to agree with alt method and keep deg of pol 1 or less    
    L_matrix(3,3) =  real_part(s_vv_s_hhs)-real_part(s_hv_s_vhs)

    
;find L-matrix vs bin

;note on notation: s_hhs=s_hh*
     
    s_vv2 = abs(c_matrix_vs_bin(0,0,*))
    s_vv_s_vhs =c_matrix_vs_bin(1,0,*)
    s_vv_s_hvs =c_matrix_vs_bin(2,0,*)
    s_vv_s_hhs = c_matrix_vs_bin(3,0,*)

    s_vh_s_vvs = c_matrix_vs_bin(0,1,*)
    s_vh2 =  abs(c_matrix_vs_bin(1,1,*))
    s_vh_s_hvs = c_matrix_vs_bin(2,1,*)
    s_vh_s_hhs = c_matrix_vs_bin(3,1,*)

    s_hv_s_vvs = c_matrix_vs_bin(0,2,*)
    s_hv_s_vhs =c_matrix_vs_bin(1,2,*)    
    s_hv2 = abs(c_matrix_vs_bin(2,2,*))
    s_hv_s_hhs =c_matrix_vs_bin(3,2,*)
    
    s_hh_s_vvs = c_matrix_vs_bin(0,3,*)
    s_hh_s_vhs =c_matrix_vs_bin(1,3,*)    
    s_hh_s_hvs = c_matrix_vs_bin(2,3,*)
    s_hh2 = abs(c_matrix_vs_bin(3,3,*))        

    L_matrix_vs_bin(0,0,*) = .5*(s_vv2+s_hh2+s_hv2+s_vh2)
    L_matrix_vs_bin(1,0,*) = .5*(s_vv2-s_hh2+s_hv2-s_vh2)
    L_matrix_vs_bin(2,0,*) = real_part(s_hv_s_hhs)+real_part(s_vv_s_vhs)
    L_matrix_vs_bin(3,0,*) = imaginary(s_hh_s_hvs)+imaginary(s_vh_s_vvs)
    
    L_matrix_vs_bin(0,1,*) = .5*(s_vv2-s_hh2-s_hv2+s_vh2)
    L_matrix_vs_bin(1,1,*) = .5*(s_vv2+s_hh2-s_hv2-s_vh2)
    L_matrix_vs_bin(2,1,*) = real_part(s_vv_s_vhs)-real_part(s_hh_s_hvs)
    L_matrix_vs_bin(3,1,*) = imaginary(s_vh_s_vvs)+imaginary(s_hv_s_hhs)
    
    L_matrix_vs_bin(0,2,*) = real_part(s_vh_s_hhs)+real_part(s_vv_s_hvs)
    L_matrix_vs_bin(1,2,*) = real_part(s_vv_s_hvs)-real_part(s_hh_s_vhs)
    L_matrix_vs_bin(2,2,*) = real_part(s_vv_s_hhs)+real_part(s_vh_s_hvs)
    L_matrix_vs_bin(3,2,*) = imaginary(s_hh_s_vvs)+imaginary(s_vh_s_hvs)    
    
    L_matrix_vs_bin(0,3,*) = imaginary(s_vh_s_hhs)+imaginary(s_vv_s_hvs)
    L_matrix_vs_bin(1,3,*) = imaginary(s_vv_s_hvs)-imaginary(s_vh_s_hhs)
    L_matrix_vs_bin(2,3,*) = -(imaginary(s_hh_s_vvs)+imaginary(s_hv_s_vhs));neg. sign added 8/1/19 to keep degree of pol <1 and agree with alt method
    L_matrix_vs_bin(3,3,*) =  real_part(s_vv_s_hhs)-real_part(s_hv_s_vhs)

 end                            ;mueller_compute

PRO mueller_compute_alt,smat,n_bins,L_matrix,L_matrix_vs_bin,n_blocks_per_group
;smat = complexarr(2,2,n_blocks_per_group,n_bins)
  r_mat=complexarr(4,4)
  ;in IDL, r_mat(i,j): i=column; j=row
r_mat(0,0)=1
r_mat(0,1)=1
r_mat(1,0)=1
r_mat(1,1)=-1
r_mat(2,2)=1
r_mat(3,2)=1
r_mat(2,3)=complex(0,-1)
r_mat(3,3)=complex(0,1)

r_mat_inv=invert(r_mat)

w_mat=complexarr(4,4)

;find L_matrix averaged over all ranges
;in IDL, w_mat(i,j): i=column; j=row

w_mat(0,0)=total(conj(smat(0,0,*,*))*smat(0,0,*,*))/n_blocks_per_group
w_mat(1,0)=total(conj(smat(1,0,*,*))*smat(1,0,*,*))/n_blocks_per_group
w_mat(2,0)=total(conj(smat(1,0,*,*))*smat(0,0,*,*))/n_blocks_per_group
w_mat(3,0)=total(conj(smat(0,0,*,*))*smat(1,0,*,*))/n_blocks_per_group

w_mat(0,1)=total(conj(smat(0,1,*,*))*smat(0,1,*,*))/n_blocks_per_group
w_mat(1,1)=total(conj(smat(1,1,*,*))*smat(1,1,*,*))/n_blocks_per_group
w_mat(2,1)=total(conj(smat(1,1,*,*))*smat(0,1,*,*))/n_blocks_per_group
w_mat(3,1)=total(conj(smat(0,1,*,*))*smat(1,1,*,*))/n_blocks_per_group

w_mat(0,2)=total(conj(smat(0,1,*,*))*smat(0,0,*,*))/n_blocks_per_group
w_mat(1,2)=total(conj(smat(1,1,*,*))*smat(1,0,*,*))/n_blocks_per_group
w_mat(2,2)=total(conj(smat(1,1,*,*))*smat(0,0,*,*))/n_blocks_per_group
w_mat(3,2)=total(conj(smat(0,1,*,*))*smat(1,0,*,*))/n_blocks_per_group

w_mat(0,3)=total(conj(smat(0,0,*,*))*smat(0,1,*,*))/n_blocks_per_group
w_mat(1,3)=total(conj(smat(1,0,*,*))*smat(1,1,*,*))/n_blocks_per_group
w_mat(2,3)=total(conj(smat(1,0,*,*))*smat(0,1,*,*))/n_blocks_per_group
w_mat(3,3)=total(conj(smat(0,0,*,*))*smat(1,1,*,*))/n_blocks_per_group


;the ## operator multiplies rows of first matrix by columns of second matrix
L_matrix=real_part(r_mat##(W_mat##r_mat_inv))
;find L-matrix vs bin
L_matrix_vs_bin=fltarr(4,4,n_bins)
;in IDL, w_mat(i,j): i=column; j=row
for ibin=0,n_bins-1 do begin
   w_mat(0,0)=total(conj(smat(0,0,*,ibin))*smat(0,0,*,ibin))/n_blocks_per_group
   w_mat(1,0)=total(conj(smat(1,0,*,ibin))*smat(1,0,*,ibin))/n_blocks_per_group
   w_mat(2,0)=total(conj(smat(1,0,*,ibin))*smat(0,0,*,ibin))/n_blocks_per_group
   w_mat(3,0)=total(conj(smat(0,0,*,ibin))*smat(1,0,*,ibin))/n_blocks_per_group

   w_mat(0,1)=total(conj(smat(0,1,*,ibin))*smat(0,1,*,ibin))/n_blocks_per_group
   w_mat(1,1)=total(conj(smat(1,1,*,ibin))*smat(1,1,*,ibin))/n_blocks_per_group
   w_mat(2,1)=total(conj(smat(1,1,*,ibin))*smat(0,1,*,ibin))/n_blocks_per_group
   w_mat(3,1)=total(conj(smat(0,1,*,ibin))*smat(1,1,*,ibin))/n_blocks_per_group

   w_mat(0,2)=total(conj(smat(0,1,*,ibin))*smat(0,0,*,ibin))/n_blocks_per_group
   w_mat(1,2)=total(conj(smat(1,1,*,ibin))*smat(1,0,*,ibin))/n_blocks_per_group
   w_mat(2,2)=total(conj(smat(1,1,*,ibin))*smat(0,0,*,ibin))/n_blocks_per_group
   w_mat(3,2)=total(conj(smat(0,1,*,ibin))*smat(1,0,*,ibin))/n_blocks_per_group

   w_mat(0,3)=total(conj(smat(0,0,*,ibin))*smat(0,1,*,ibin))/n_blocks_per_group
   w_mat(1,3)=total(conj(smat(1,0,*,ibin))*smat(1,1,*,ibin))/n_blocks_per_group
   w_mat(2,3)=total(conj(smat(1,0,*,ibin))*smat(0,1,*,ibin))/n_blocks_per_group
   w_mat(3,3)=total(conj(smat(0,0,*,ibin))*smat(1,1,*,ibin))/n_blocks_per_group

 
   L_matrix_vs_bin(*,*,ibin)=real_part(r_mat##(W_mat##r_mat_inv))
endfor



 end                            ;mueller_compute_alt
PRO amat_fmat_compute,spec,scatvars,configvars,calvars,range_gate_spacing,n_blocks,range,ngates,ainv,finv,corner_range,$
total_corner_power_vv,total_corner_power_hh,scan_index,sweep_count,elevation,current_calibration_loop_power,base_filename,calfile

;
;compute amat and fmat from isotropic data
;see: /projects/180517-Manitoba_Ku-Ka-Scats/analysis/polarimetric_calibration_methodology.pdf
; 


window,0,re=2,xsize=450,ysize=300,xpos=240,ypos=225
window,1,re=2,xsize=450,ysize=300,xpos=800,ypos=225
window,2,re=2,xsize=450,ysize=300,xpos=240,ypos=550
window,3,re=2,xsize=450,ysize=300,xpos=800,ypos=550
gate_spacing=range(1)-range(0)
sample_range_width=0.025

if max(scan_index) ne 1 then begin
   n_blocks_ave=n_blocks ;stare mode-use only with the dual antenna Ku-Ka scat in stare mode 
   start_block=0
endif else begin
   index_scan=where(scan_index eq 1)
   max_el_for_cal=16.0;degrees
   index_low_el=where(elevation(index_scan) lt max_el_for_cal,count_low_el)
   if count_low_el eq 0 then begin
      print,'lowest elevation angle must be below 10 degrees for calibration'
      stop
   endif
   start_block=index_scan(index_low_el(0));fixed this on 10/21/19
   n_blocks_ave=count_low_el
endelse

avespec_vv=fltarr(ngates)
avespec_vh=fltarr(ngates)
avespec_hv=fltarr(ngates)
avespec_hh=fltarr(ngates) 

for iblock=start_block,start_block+n_blocks_ave-1 do begin ;fixed this on 10/21/19...",start_block+n_blocks_ave'
   avespec_vv=avespec_vv+abs(spec(iblock,0,*)^2)
   avespec_vh=avespec_vh+abs(spec(iblock,1,*)^2)
   avespec_hv=avespec_hv+abs(spec(iblock,2,*)^2)
   avespec_hh=avespec_hh+abs(spec(iblock,3,*)^2)
endfor
avespec_vv=avespec_vv/n_blocks_ave
avespec_vh=avespec_vh/n_blocks_ave
avespec_hv=avespec_hv/n_blocks_ave
avespec_hh=avespec_hh/n_blocks_ave
avespec_vh=avespec_vh/n_blocks_ave
sumspec=(avespec_hh+avespec_vv)/2.0
sumspec_xpol=(avespec_vh+avespec_hv)/2.0

index=where(range gt 0)
gate0range=index(0)
index_peak=where(sumspec(gate0range:ngates-1) eq max(sumspec(gate0range:ngates-1)))
index_peak_xpol=where(sumspec_xpol(gate0range:ngates-1) eq max(sumspec_xpol(gate0range:ngates-1)))
read,'enter 1 to use copol peak, 2 to use xpol peak for calibration: ',ipeak

if ipeak eq 1 then gate_peak=gate0range+index_peak(0) else gate_peak=gate0range+index_peak_xpol(0) 
testpeak =range(gate_peak)-sample_range_width/2.

new_search:
index=where(range gt testpeak)
gate0=index(0)
n_bins=ceil(sample_range_width/range_gate_spacing)+1
if n_bins gt ngates-gate0 then n_bins=ngates-gate0
gate1=gate0+n_bins-1
gatemaxplot=gate1+40
if gatemaxplot gt ngates-1 then gatemaxplot=ngates-1
wset,2
plot,range(gate0range:gatemaxplot),10*alog10(sumspec(gate0range:gatemaxplot)),xs=1,xtitle='range (m)',ytitle='power (dB)',title='cal taken in red highlighted region',charsize=1.5

oplot,range(gate0range:gatemaxplot),10*alog10(sumspec_xpol(gate0range:gatemaxplot)),color=80
;oplot,range(gate0range:gatemaxplot),10*alog10(avespec_hh(gate0range:gatemaxplot)),color=25
oplot,range(gate0:gate1),10*alog10(sumspec(gate0:gate1)),color=160,thick=2

vmat=complexarr(2,2,n_blocks_ave,n_bins)
 for iblock=start_block,start_block+n_blocks_ave-1 do begin ;fixed this on 10/21/19...",start_block+n_blocks_ave'
    iblock0=iblock-start_block
  vmat(0,0,iblock0,*) = spec(iblock,0,gate0:gate1) ;vmat(0,0) is the V V return
  vmat(1,0,iblock0,*) = spec(iblock,2,gate0:gate1)       ;vmat(1,0) is the V pol return when transmitting H-changed 7/13/19 from vmat(0,1)
  vmat(0,1,iblock0,*) = spec(iblock,1,gate0:gate1)       ;vmat(0,1) is the H pol return when transmitting V-changed 7/13/19 from vmat(1,0)
  vmat(1,1,iblock0,*) = spec(iblock,3,gate0:gate1)       ;vmat(1,1) is the H H return
endfor

covariance_vvhh=complexarr(n_bins)
covariance_vhhv=complexarr(n_bins)
vv_over_hh_power=fltarr(n_bins)
vh_over_hv_power=fltarr(n_bins)
for ibin=0,n_bins-1 do begin
   covariance_vvhh(ibin)=mean(vmat(0,0,*,ibin)*conj(vmat(1,1,*,ibin)))
   covariance_vhhv(ibin)=mean(vmat(0,1,*,ibin)*conj(vmat(1,0,*,ibin)))
   vv_over_hh_power(ibin)=mean(abs(vmat(0,0,*,ibin)^2))/mean(abs(vmat(1,1,*,ibin)^2))
   vh_over_hv_power(ibin)=mean(abs(vmat(0,1,*,ibin)^2))/mean(abs(vmat(1,0,*,ibin)^2))
endfor

thetavec=atan(covariance_vvhh,/phase)
phivec=atan(covariance_vhhv,/phase)
wset,0
plot,range(gate0:gate1),thetavec/!dtor,xtitle='range (m)',ytitle='phase (deg)',title='VV/HH phase',charsize=1.5,yrange=[-180,180]
wset,1
plot,range(gate0:gate1),phivec/!dtor,xtitle='range (m)',ytitle='phase (deg)',title='VH/HV phase',charsize=1.5,yrange=[-180,180]
wset,3
plot,range(gate0:gate1),10*alog10(VV_over_HH_power),xtitle='range (m)',ytitle='dB',title='VV/HH power ratio',charsize=1.5,yrange=[-6,6]
theta=mean(thetavec)
phi=mean(phivec)
VV_over_HH_pwr_mean=mean(VV_over_hh_power)
VH_over_HV_pwr_mean=mean(vh_over_hv_power)
print
print,'mean angle vv/hh (deg): ',theta/!dtor
print,'mean angle vh/hv (deg): ',phi/!dtor
print,'mean VV over HH power (dB): ',10*alog10(VV_over_HH_pwr_mean)

goto,jump_changes
   print,'enter "l" to move peak search range to left, "r" to move right, or "s" to enter peak range, otherwise any key to accept range (q to quit), m to set new sample range width:'
   dum = get_kbrd(1)
   IF dum EQ 'q' THEN stop
   IF dum EQ 'Q' THEN stop
   testpeak0=testpeak
;testpeak is the guess of peak location in meters
   IF dum EQ 'l' THEN BEGIN 
      testpeak = testpeak-sample_range_width/4   ;meters;
      GOTO,new_search
   ENDIF 
   IF dum EQ 'r' THEN BEGIN 
      testpeak = testpeak+sample_range_width/4
      GOTO,new_search
   ENDIF 
   IF dum EQ 's' THEN BEGIN 
      read,'enter range of peak: ',testpeak
      GOTO,new_search
   ENDIF
   IF dum EQ 'm' THEN BEGIN
      print,'current sample range width (m): ',sample_range_width
      print,'minimum range width (m): ',gate_spacing
      read,'enter sample range width: ',sample_range_width
      if sample_range_width lt gate_spacing then begin
         print,'sample range width set to gate spacing'
         sample_range_width=gate_spacing
         waitkey,dum
      endif

      GOTO,new_search
  ENDIF 

jump_changes:
mag_beta=(VV_over_HH_pwr_mean/VH_over_HV_pwr_mean)^.25
mag_alpha=(VV_over_HH_pwr_mean*VH_over_HV_pwr_mean)^.25
phi_alpha=(theta-phi)/2.
phi_beta=(theta+phi)/2.

alpha=mag_alpha*exp(complex(0.0,phi_alpha))
beta=mag_beta*exp(complex(0.0,phi_beta))

ainv=complexarr(2,2)
finv=complexarr(2,2)

ainv(0,0)=1.0
ainv(1,1)=alpha

finv(0,0)=1.0
finv(1,1)=beta

amat=invert(ainv)
fmat=invert(finv)
print
print,'ainv: '
print,ainv
print
print,'finv: '
print,finv
print
irepeat=1
;read,'enter 1 to accept 0 to change sample range: ',irepeat
if irepeat eq 0 then goto,new_search
calfile=base_filename+'.cal'
print
print
print
print,'creating new calibration file, name: ',calfile
print
print
;read,'input new calibration file name: ',calfile

total_corner_power_vv=10^(calvars.corner_reflector_vv_power_dbm/10.)
total_corner_power_hh=10^(calvars.corner_reflector_hh_power_dbm/10.)
corner_range=calvars.corner_reflector_range_m

print,'total_corner_power_vv (dBm): ',calvars.corner_reflector_vv_power_dbm
print,'total_corner_power_hh (dBm): ',calvars.corner_reflector_hh_power_dbm
print,'corner range (m): ',corner_range
  openw,1,calfile
  printf,1,ainv,finv,current_calibration_loop_power,corner_range,total_corner_power_vv,total_corner_power_hh,amat,fmat
  close,1
waitkey,dum
end;amat_fmat_compute

PRO equalize,scatvars,padfac,H_k_norm
                                ;this subroutine computes the
                                ;equalization function required to
                                ;compensate for gain loss as a
                                ;function of Doppler frequency due to
                                ;block averaging.  See
                                ;/home/mead/Manitoba_Ku-Ka_Scats/idl/block_averaging_freq_response.pro
                                ;or online discussions of frequency
                                ;response of running average filter
  if scatvars.decimation_mode eq 0 then begin 
     Lfac=long(scatvars.decimation)
     N_predec=long(scatvars.n_gates*Lfac) ;predecimation number of samples
     N_fft=long(scatvars.n_gates*padfac)
     delta_k=1e-10;added 12/19/19+next line modified
     kvec=lindgen(N_predec)+delta_k
     H_k=1.0/Lfac*(1-exp(dcomplex(0.0,-2*!pi*kvec*Lfac/N_predec)))/(1-exp(complex(0.0,-2*!pi*kvec/N_predec)))
     H_k_norm0=abs(H_k)/max(abs(H_k))
     H_k_norm1=H_k_norm0(0:N_predec/(Lfac)-1)
     h_k_norm=congrid(H_k_norm1,N_fft)
  endif else begin
     h_k_norm=1.0
  endelse

 ; stop
end

PRO process_pol_data_in_footprint,configvars,c_matrix_vs_bin,L_matrix_vs_bin,range,gate0,gate1,n_bins
      cor_coe_vs_bin = c_matrix_vs_bin(0,3,*)/sqrt(c_matrix_vs_bin(0,0,*)*c_matrix_vs_bin(3,3,*))   
       cor_coe_xpol_vs_bin=c_matrix_vs_bin(1,2,*)/sqrt(c_matrix_vs_bin(1,1,*)*c_matrix_vs_bin(2,2,*))      
       
       rho_hv_vs_bin = abs(cor_coe_vs_bin)
       rho_hv_xpol_vs_bin = abs(cor_coe_xpol_vs_bin)
       
       phase_hv_vs_bin = atan(imaginary(cor_coe_vs_bin),real_part(cor_coe_vs_bin))/!dtor
       phase_hv_xpol_vs_bin = atan(imaginary(cor_coe_xpol_vs_bin),real_part(cor_coe_xpol_vs_bin))/!dtor
       
       rel_mag_vs_bin=(c_matrix_vs_bin(0,0,*)+c_matrix_vs_bin(3,3,*))/max((c_matrix_vs_bin(0,0,*)+c_matrix_vs_bin(3,3,*)))
       rel_mag_xpol_vs_bin=(c_matrix_vs_bin(1,1,*)+c_matrix_vs_bin(2,2,*))/max((c_matrix_vs_bin(0,0,*)+c_matrix_vs_bin(3,3,*)))
       LDR_vs_bin=(c_matrix_vs_bin(1,1,*)+c_matrix_vs_bin(2,2,*))/(c_matrix_vs_bin(0,0,*)+c_matrix_vs_bin(3,3,*))
;plot copol data vs bin
       
       wset,2
       plot,range(gate0:gate1),10*alog10(rel_mag_vs_bin),xtitle='range (m)',title='relative magnitude (copol: blk; xpol red) ',ytitle='dB',charsize=1.4,xs=1,yr=[-60,0],ys=1
       oplot,range(gate0:gate1),10*alog10(rel_mag_xpol_vs_bin),color=160
       wset,0
       plot,range(gate0:gate1),rho_hv_vs_bin,xtitle='range (m)',ytitle='magnitude',title='mag correlation coe. (copol blk; xpol red)',yr=[0,1.1],ys=1,charsize=1.4,xs=1
       oplot,range(gate0:gate1),rho_hv_xpol_vs_bin,color=160
       xvec=[range(gate0),range(gate1)]
       yvec=[1,1]
       oplot,xvec,yvec,line=2
       wset,1
       plot,range(gate0:gate1),phase_hv_vs_bin,xtitle='range (m)',ytitle='phase (deg)',title='phase correlation coefficient',yr=[-180,180],charsize=1.4,ys=1,xs=1
       wset,3
       plot,float(cor_coe_vs_bin),imaginary(cor_coe_vs_bin),psym=2,yrange=[-1,1],xr=[-1,1],title='complex correlation coefficient',charsize=1.4
       oplot,float(cor_coe_vs_bin(0:n_bins/20)),imaginary(cor_coe_vs_bin(0:n_bins/20)),psym=2,color=160
       oplot,float(cor_coe_vs_bin(.95*n_bins:n_bins-1)),imaginary(cor_coe_vs_bin(.95*n_bins:n_bins-1)),psym=2,color=60
;
;plot additional xpol data
;
       waitkey,dum
       wset,1
       plot,range(gate0:gate1),phase_hv_xpol_vs_bin ,xtitle='range (m)',ytitle='phase (deg)',title='phase xpol correlation coefficient',yr=[-180,180],charsize=1.4,ys=1,xs=1
       wset,3
       plot,float(cor_coe_xpol_vs_bin),imaginary(cor_coe_xpol_vs_bin),psym=2,yrange=[-1,1],xr=[-1,1],title='xpol complex correlation coefficient',charsize=1.4
       oplot,float(cor_coe_xpol_vs_bin(0:n_bins/20)),imaginary(cor_coe_xpol_vs_bin(0:n_bins/20)),psym=2,color=160
       oplot,float(cor_coe_xpol_vs_bin(.95*n_bins:n_bins-1)),imaginary(cor_coe_xpol_vs_bin(.95*n_bins:n_bins-1)),psym=2,color=60
       
       wset,2
       plot,range(gate0:gate1),10*alog10(LDR_vs_bin),xtitle='range (m) ',title='Linear Depolarization Ratio ',ytitle='LDR (dB)',charsize=1.4,xs=1,yr=[-40,2],ys=1
       xvec=[range(gate0),range(gate1)]
       
       yvec=[0,0]
       oplot,xvec,yvec,line=2
       wset,3
       plot,LDR_vs_bin,cor_coe_vs_bin,psym=2,xr=[0,1.2],yr=[0,1],xtitle='LDR',ytitle='rho_hv',title='rhohv vs LDR monostatic/bistatic model',charsize=1.5
       xvec=[0,.5]
       yvec=[1,0]
       oplot,xvec,yvec,line=2
       xvec=[0,1]
       oplot,xvec,yvec,line=2
       waitkey,dum
       proc_thresh_min_dB=min([configvars.proc_thresh_left_dB,configvars.proc_thresh_right_dB])
       if configvars.i_pol_signature eq 1 then begin
          i_multi_bin=1
          control_plots=1
          iline_dum=0
          for ibin=0,n_bins-1,configvars.smoothfac do begin
             wset,3
             plot,range(gate0:gate1),10*alog10(rel_mag_vs_bin),xtitle='range (m)',title='co-pol rel. mag (blk); LDR (green): ',ytitle='dB',charsize=1.4,xs=1,yr=[proc_thresh_min_dB-10,0],ys=1
             oplot,range(gate0:gate1),10*alog10(LDR_vs_bin),color=80
             xvec=[range(gate0+ibin),range(gate0+ibin)]
             yvec=[proc_thresh_min_dB-10,0]
             oplot,xvec,yvec,color=160
             L_matrix_bin=L_matrix_vs_bin(*,*,ibin)
waitkey,dum
             polsig,L_matrix_bin,iline_dum,control_plots,ibin,i_multi_bin
          endfor                ;ibin
       endif                    ;i_pol_signature
end;process_pol_data_in_footprint

PRO process_polarimetric_data,scatvars,configvars,calvars,range,spec,n_blocks_per_group,ngates,smoothfac,c_matrix_vec,L_matrix_vec,rho_hv_vec,phase_hv_deg_vec,rdepol_db_vec,total_power,range_peak_signal,range_centroid_signal,igroup,group_index,gate_offset,gate_plot_max,gate_peak,gate_peak_cal,ainv,finv,scan_index,line_el,peak_power
  
c_matrix=complexarr(4,4)
L_matrix=fltarr(4,4)

;set avespec to zero each call of subroutine
 avespec_vv = fltarr(ngates)  ;average power spectrum (power spectrum=range profile of signal power)
 avespec_hv =  fltarr(ngates)
 avespec_vh =  fltarr(ngates)
 avespec_hh =  fltarr(ngates)
 avespec_cal=fltarr(ngates)
;average range profiles 
  print
  print,'processing group/line number: ',igroup
  print
  for jsweep=0,n_blocks_per_group-1 do begin
      isweep = group_index(igroup)+jsweep
      avespec_vv=avespec_vv+abs(spec(isweep,0,*)^2)
      avespec_vh=avespec_vh+abs(spec(isweep,1,*)^2)
      avespec_hv=avespec_hv+abs(spec(isweep,2,*)^2)
      avespec_hh=avespec_hh+abs(spec(isweep,3,*)^2)
      avespec_cal=avespec_cal+abs(spec(isweep,4,*)^2)
  endfor                       ;jsweep
   avespec_vv = avespec_vv/n_blocks_per_group
   avespec_hv = avespec_hv/n_blocks_per_group
   avespec_vh = avespec_vh/n_blocks_per_group
   avespec_hh = avespec_hh/n_blocks_per_group
   avespec_cal = avespec_cal/n_blocks_per_group
 
   sumspec=(avespec_vv+avespec_hh+avespec_vh+avespec_hv)/4.
;
;n_bins=number of bins (spaced by range_gate_spacing for which a covariance matrix is
;computed 
;
testpeak =scatvars.pedestal_height_m
index=where(range gt testpeak,count)
search_margin=ngates/50
index_start_search=index(0)-search_margin
if index_start_search lt ngates/30 then index_start_search=ngates/30;avoid DC peak

proc_thresh_right=10^(configvars.proc_thresh_right_dB/10)
proc_thresh_left=10^(configvars.proc_thresh_left_dB/10)
index=where(sumspec(index_start_search:ngates-1) eq max(sumspec(index_start_search:ngates-1)))
gate_peak=index(0)+index_start_search


testpeak_cal=-1.0
index=where(range gt testpeak_cal,count)
index_start_search_cal=index(0)
index_cal=where(avespec_cal(index_start_search_cal:ngates-1) eq max(avespec_cal(index_start_search_cal:ngates-1)))
gate_peak_cal=index_cal(0)+index_start_search_cal


iloop_left=0
loop0=0
pre_peak_gate=fix(gate_peak*.8);recoded finding of peak region on 4/27/2020
retry_left:
index_lo_left=where(sumspec(pre_peak_gate:gate_peak)/sumspec(gate_peak) lt proc_thresh_left,count)
print,'range(index_lo_left(0)): ',range(index_lo_left(0)+pre_peak_gate)
if count eq 0 then begin 
   proc_thresh_left = 2*proc_thresh_left
   if proc_thresh_left gt 1 or loop0 gt 0 then proc_thresh_left=0.5+loop0/10
   loop0=loop0+1
   if iloop_left gt 10 then begin
      print,'cannot find signal region left of peak'
      stop
   endif
   iloop_left=iloop_left+1
   goto,retry_left
endif
gate0=pre_peak_gate+max(index_lo_left)

if gate0 gt ngates-2 then gate0=ngates-2
if gate0 lt 0 then gate0=0

iloop_right=0
retry_right:
index_lo_right=where(sumspec(gate_peak:ngates-1)/max(sumspec(gate_peak:ngates-1)) lt proc_thresh_right,count)
if count eq 0 then begin 
   proc_thresh_right = 10*proc_thresh_right
   if proc_thresh_right gt 1 then proc_thresh_right=0.5
   if iloop_right gt 10 then begin
      print,'cannot find signal region right of peak'
      stop
   endif
   iloop_right=iloop_right+1 
   goto,retry_right
endif
gate1=gate_peak+index_lo_right(0)
if gate1-gate0 lt 1 then gate1=gate0+1

if configvars.i_corner_process eq 1 then begin
   gate0=gate_peak-calvars.n_summed_gates/2
   gate1=gate0+calvars.n_summed_gates-1
endif

n_bins=gate1-gate0+1

if n_bins le smoothfac then begin; le changed from lt on 4/13/2020
   print,'n_bins: ',n_bins
   print,'n_bins less than smoothfac'
   print,'reducing smoothing factor to n_bins/2'
   smoothfac=n_bins/2
endif


range_peak_signal(igroup) =range(gate_peak)
range_centroid_vv = total(range(gate0:gate1)*avespec_vv(gate0:gate1))/total(avespec_vv(gate0:gate1))
range_centroid_hh = total(range(gate0:gate1)*avespec_hh(gate0:gate1))/total(avespec_hh(gate0:gate1)) 
range_centroid_signal(igroup) = (range_centroid_vv+range_centroid_hh)/2.

c_matrix_vs_bin = complexarr(4,4,n_bins) ;covariance matrix at each range within beam footprint
L_matrix_vs_bin = fltarr(4,4,n_bins) ;covariance matrix at each range within beam footprint
vmat = complexarr(2,2,n_blocks_per_group,n_bins)
smat = complexarr(2,2,n_blocks_per_group,n_bins)
 
total_power_vv = total(avespec_vv(gate0:gate1))
total_power_hv = total(avespec_hv(gate0:gate1))
total_power_vh = total(avespec_vh(gate0:gate1))
total_power_hh = total(avespec_hh(gate0:gate1))

total_power(igroup,0)=total_power_vv
total_power(igroup,1)=total_power_hv
total_power(igroup,2)=total_power_vh
total_power(igroup,3)=total_power_hh

peak_power_vv=max(avespec_vv(gate0:gate1))
peak_power_hv=max(avespec_hv(gate0:gate1))
peak_power_vh=max(avespec_vh(gate0:gate1))
peak_power_hh=max(avespec_hh(gate0:gate1))

peak_power(igroup,0)=peak_power_vv
peak_power(igroup,1)=peak_power_hv
peak_power(igroup,2)=peak_power_vh
peak_power(igroup,3)=peak_power_hh

print,'total_power_vv (dB): ',10*alog10(total_power_vv)
print,'total_power_hh (dB): ',10*alog10(total_power_hh)
print,'total_power_vh (dB): ',10*alog10(total_power_vh)
print,'total_power_hv (dB): ',10*alog10(total_power_hv)

   wset,0
   erase
   wset,1
   erase
   wset,3
   erase
wset,1
yrmin=10*alog10(min(avespec_vv(gate_offset:ngates-1)))-3
yrmax=10*alog10(max([avespec_vv(gate_offset:ngates-1),avespec_hh(gate_offset:ngates-1)]))+3
gate_plot=fix(gate1*1.05)
if gate_plot gt ngates then gate_plot=float(.2*ngates)
plot,range(gate_offset:gate_plot_max),10*alog10(avespec_vv(gate_offset:gate_plot_max)),xtitle='range (m)',ytitle='relative power (dB)',title='ave. power VV blk; HH red; VH blue; HV grn',charsize=1.4,yrange=[yrmin,yrmax],xs=1,ys=1
oplot,range(gate_offset:gate_plot_max),10*alog10(avespec_hh(gate_offset:gate_plot_max)),color=160
oplot,range(gate_offset:gate_plot_max),10*alog10(avespec_vh(gate_offset:gate_plot_max)),color=60 
oplot,range(gate_offset:gate_plot_max),10*alog10(avespec_hv(gate_offset:gate_plot_max)),color=40
db_str_vv=strmid(strtrim((fix(1000*alog10(peak_power_vv)))/100.,2),0,5)+' dB'
db_str_hh=strmid(strtrim((fix(1000*alog10(peak_power_hh)))/100.,2),0,5)+' dB'
xyouts,range(gate_plot),yrmin+.9*(yrmax-yrmin),'VV peak: '+db_str_vv,charsize=1.1,alignment=0
xyouts,range(gate_plot),yrmin+.82*(yrmax-yrmin),'HH peak: '+db_str_hh,charsize=1.1,alignment=0
if max(scan_index) eq 1 then begin
   el_ang_str=strtrim(round(line_el),2)
   xyouts,range(gate_plot),yrmin+.74*(yrmax-yrmin),'elevation angle (deg): '+el_ang_str,charsize=1.1,alignment=0.0
endif
rvec1 = [range(gate0),range(gate0)]             
rvec2 = [range(gate1),range(gate1)]             
yvec = [yrmin,yrmax]
oplot,rvec1,yvec,color=1000
oplot,rvec2,yvec,color=1000

wset,0
yrmin=max(10*alog10(avespec_vv(gate0:gate1)))+min([configvars.proc_thresh_left_dB,configvars.proc_thresh_right_dB])-5
yrmax=max(10*alog10(avespec_vv(gate0:gate1)))+5
plot,range(gate0:gate1),10*alog10(avespec_vv(gate0:gate1)),xtitle='range (m)',ytitle='relative power (dB)',title='ave. power VV blk; HH red; VH blue; HV grn',charsize=1.4,xs=1,ys=1,yr=[yrmin,yrmax]
oplot,range(gate0:gate1),10*alog10(avespec_hh(gate0:gate1)),color=160
oplot,range(gate0:gate1),10*alog10(avespec_vh(gate0:gate1)),color=60 
oplot,range(gate0:gate1),10*alog10(avespec_hv(gate0:gate1)),color=4
;oplot,range(gate0:gate1),10*alog10(avespec_hh(gate0:gate1)*range(gate0:gate1)^2/range(gate_peak)^2),color=25,thick=2,line=2
gate_plot1=.4*fix(gate0+gate1)
xyouts,mean(range(gate0:gate1)),yrmin+.14*(yrmax-yrmin),'data selected for processing ',charsize=1.5,alignment=.5
waitkey,dum
;compute scattering matrix for target under observation; first compute
;voltage matrix (vmat), then calibrate to get scattering matrix (smat)
;
;vmat is temporary-reused for each igroup
;
;note: v_mat(i,j) is the ith column, jth row
   for jsweep=0,n_blocks_per_group-1 do begin
       isweep = group_index(igroup)+jsweep
        vmat(0,0,jsweep,*) = spec(isweep,0,gate0:gate1) ;vmat(0,0) is the V V return
        vmat(1,0,jsweep,*) = spec(isweep,2,gate0:gate1) ;vmat(1,0) is the V pol return when transmitting H-changed 7/13/19 from vmat(0,1)
        vmat(0,1,jsweep,*) = spec(isweep,1,gate0:gate1) ;vmat(0,1) is the H pol return when transmitting V-changed 7/13/19 from vmat(1,0)
        vmat(1,1,jsweep,*) = spec(isweep,3,gate0:gate1) ;vmat(1,1) is the H H return
   endfor
;compute complex scattering matrix, Smat; data becomes calibrated by applying ainv and finv
    FOR j=0,n_blocks_per_group-1 DO BEGIN
        FOR k=0,n_bins-1 DO BEGIN 
             Smat(*,*,j,k) = ainv##vmat(*,*,j,k)##finv;this step computes calibrated 2x2 scattering matrix at each range bin within gate0:gate1       
         ENDFOR
    ENDFOR  

;compute covariance of S_vv and S_hh
    npoints = n_blocks_per_group*n_bins
;
;display polarization ellipse scatter plot
;
   wset,1

   if configvars.i_pol_scat eq 1 then   polscatter,smat,npoints,line_el

;reform 2x2 smat to 4x1 vector form at each sweep and bin
    smatvec = reform(smat,4,n_blocks_per_group,n_bins)

;compute 4x4 covariance matrix containing calibrated covariances of
;all combinations of scattering matrix elements for each line

    FOR j=0,3 DO BEGIN
        FOR k=0,3 DO BEGIN
                                ;smatvec(4 elements,time,distance)
                                ;reversed indices in next line on 7/13/19 to comply with column/row indexing
            c_matrix(j,k) = total(smatvec(k,*,*)*conj(smatvec(j,*,*)))/n_blocks_per_group;don't divide by n_bins here...we are summing over all range bins; divide by group_index since we're averaging power in time 
         ENDFOR 
     ENDFOR
for m=0,n_bins-1 do begin;range gate dimension
    FOR j=0,3 DO BEGIN;dimension of 4x4 covariance matrix
        FOR k=0,3 DO BEGIN;dimension of 4x4 covariance matrix
           ;averaging over time
            c_matrix_vs_bin(j,k,m) = total(smatvec(k,*,m)*conj(smatvec(j,*,m)))/n_blocks_per_group
         ENDFOR                 ;k 
     ENDFOR;j
 endfor;m
;print,c_matrix_vs_bin(*,*,0,0)
;waitkey,dum
smooth_vec=intarr(3)
smooth_vec(*)=1
smooth_vec(2)=smoothfac
c_matrix_vs_bin=smooth(c_matrix_vs_bin,smooth_vec,/edge_truncate);smoothing over bin (range bin or gate) dimension

mueller_compute,c_matrix,c_matrix_vs_bin,L_matrix,L_matrix_vs_bin

c_matrix_vec(*,*,igroup)=c_matrix
L_matrix_vec(*,*,igroup)=L_matrix

print,'Mueller matrix averaged over all range bins: '
print,'L_matrix: '
print,L_matrix/max(L_matrix)
print
print,'L_matrix vs bin for first bin: '
print,L_matrix_vs_bin(*,*,0)/max(L_matrix_vs_bin(*,*,0))

;the following code was used to check validity of mueller_compute
;using c_matrix formulation.  These agree when smoothfac=1
goto,skip_alt
mueller_compute_alt,smat,n_bins,L_matrix_alt,L_matrix_vs_bin_alt,n_blocks_per_group
print
print,'Mueller matrix averaged over all range bins (alternate form): '
print,L_matrix_alt/max(L_matrix_alt)
print
print,'L_matrix vs bin (alternate form): '
print,L_matrix_vs_bin_alt(*,*,0)/max(L_matrix_vs_bin_alt(*,*,0))
L_matrix=L_matrix_alt
L_matrix_vs_bin=L_matrix_vs_bin_alt
skip_alt:

cor_coe = c_matrix(0,3)/sqrt(c_matrix(0,0)*c_matrix(3,3))   
cor_coe_xpol=c_matrix(1,2)/sqrt(c_matrix(1,1)*c_matrix(2,2))   
rho_hv = abs(cor_coe)
phase_hv = atan(imaginary(cor_coe),real_part(cor_coe))
rdepol = abs(c_matrix(2,1))/abs(sqrt(c_matrix(0,0)*c_matrix(3,3)))
rco = float(c_matrix(0,0)/c_matrix(3,3))
print
print,'mag co-polarized correlation coefficient of Smat: ',rho_hv
print,'phase co-polarized correlation coefficient of Smat (degrees): ',phase_hv/!dtor
print,'depolarization ratio (dB): ',10*alog10(rdepol)
print,'co-polarized ratio (rco) (dB): ',10*alog10(rco)
rho_hv_vec(igroup)=rho_hv
phase_hv_deg_vec(igroup)=phase_hv/!dtor
rdepol_dB_vec(igroup)=10*alog10(rdepol)   
wset,0
;
;display polarization signature of target averaged over all bins
;

if max(scan_index) eq 1 then begin
   control_plots=1
   i_multi_bin=0
   iline_dum=igroup
endif else begin
   control_plots=1
   line_el=0                    ;
   i_multi_bin=0
   iline_dum=0
endelse
if configvars.i_pol_signature eq 1 then begin
   polsig,L_matrix,iline_dum,control_plots,line_el,i_multi_bin
endif

for i=0,3 do begin
   wshow,i
endfor
if configvars.i_pol_data_in_footprint eq 1 then process_pol_data_in_footprint,configvars,c_matrix_vs_bin,L_matrix_vs_bin,range,gate0,gate1,n_bins

end                             ;process_polarimetric_data

PRO earth_radius,lat,R_e

;formula from Wikipedia
;Distances produced by Haversine formula now agree very closely with
;Google Earth (within a cm at 760 meters).

a_equitorial=6378137.0          ;Earth's radius at the equator
b_polar=6356752.3               ;Earth's radius at the poles

num=(a_equitorial^2*cos(lat))^2+(b_polar^2*sin(lat))^2
denom=(a_equitorial*cos(lat))^2+(b_polar*sin(lat))^2

R_e=sqrt(num/denom)

;print,'radius at latitude (deg): ',R_e,lat/!dtor

end                             ;earth_radius
PRO compute_distance,gps_latitude,gps_longitude,elapsed_time,scatvars,configvars,distance,independent_sample_index
;
;Use Haversine formula to compute distance and bearing between two lat/lon points
;reference: http://www.movable-type.co.uk/scripts/latlong.html

  lat=gps_latitude
  lon=gps_longitude
  time=elapsed_time
  nvals=n_elements(gps_latitude)
  xoff=fltarr(nvals)
  yoff=fltarr(nvals)
  distance=fltarr(nvals)
  bearing=fltarr(nvals)
  
  lat1 = double(lat(0)*!dtor)
  lon1 = double(lon(0)*!dtor)
  earth_radius,lat1,R_e
  for i=0,nvals-1 do begin
     lat2 = double(lat(i)*!dtor)
     lon2 = double(lon(i)*!dtor)
     delta_lat = (lat2-lat1)
     delta_lon = (lon2-lon1)
     
     afac = sin(delta_lat/2)^2+cos(lat1)*cos(lat2)*sin(delta_lon/2.)^2
     
     cfac = 2.0*atan(sqrt(afac),sqrt(1.-afac))
     
     distance(i) = float(R_e*cfac)
     
 ;    print,'distance between lat/long in m: ',distance(i)
     
;find bearing in radians (bearing is angular position of point from
;origin relative to true north) 
     
     bearing(i) = float(atan(sin(delta_lon)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*(cos(delta_lon))))
     
     IF bearing(i) LT 0 THEN bearing(i) = bearing(i)+2*!pi ;
     
;     print,'bearing in degrees: ',bearing(i)/!dtor   
     xoff(i)=distance(i)*sin(bearing(i))
     yoff(i)=distance(i)*cos(bearing(i))
  endfor
;smooth out distance to account for GPS short-term errors
  smoothfac=20
  if smoothfac gt nvals then smoothfac=nvals/6+1
  distance=smooth(distance,smoothfac,/edge_truncate)
  independent_sample_index0=lonarr(nvals)
  distance_ref=distance(0)
;compute velocity
  vec0=indgen(nvals-1)
  vec1=vec0+1
  velocity0=(distance(vec1)-distance(vec0))/(elapsed_time(vec1)-elapsed_time(vec0))
  velocity=[velocity0,velocity0(nvals-2)]

;compute distance to independence for finding indices of data which
;are statistically independent 

  Ku_ant_diameter=0.106;was incorrectly set to .15   
  Ka_ant_diameter=0.064;was incorrectly set to 0.09
  distance_to_independence=Ku_ant_diameter/2.0
  if scatvars.chirp_center_frequency_hz gt 20e09 then distance_to_independence=Ka_ant_diameter/2.0
;test for independent samples
  min_velocity=0.5; test for min velocity to avoid drifting GPS location from appearing as motion
  if configvars.i_proc_ind eq 0 then begin
     min_velocity =-999.
     distance_to_independence=-999.
  endif
  count_ind=0
  for i=0,nvals-1 do begin
     if abs(distance(i)-distance_ref) gt distance_to_independence and velocity(i) gt min_velocity then begin
        independent_sample_index0(count_ind)=i
        count_ind=count_ind+1
        distance_ref=distance(i)
     endif
  endfor
  independent_sample_index=independent_sample_index0(0:count_ind-1) 
  
  wset,3
  erase
  wset,0
  xmax=max(abs(xoff))
  ymax=max(abs(yoff))
  cart_max=max([xmax,ymax])*1.05
  plot,xoff,yoff,charsize=1.5,xtitle='east/west distance (m)',ytitle='north/south distance (m)',title='displacement from origin',xr=[-cart_max*2,cart_max*2],yr=[-cart_max,cart_max],ys=1,xs=1
  oplot,xoff,yoff,color=160
  wset,1
  plot,elapsed_time/60,velocity,charsize=1.5,xtitle='elapsed time (m)',ytitle='velocity (m/s)',title='velocity estimated from position derivative'
  xvec=[0,max(elapsed_time/60)]
  yvec=[min_velocity,min_velocity]
  oplot,xvec,yvec,line=2,color=160
  xyouts,mean(xvec),yvec*1.07,'velocity threshold',color=160,alignment=.5,charsize=1.5
  wset,2
  titlestr='location of independent samples shown in red'
  if configvars.i_proc_ind eq 0 then titlestr='location of all samples'
  plot,elapsed_time/60,distance,charsize=1.5,xtitle='elapsed time (m)',ytitle='along track distance (m)',title=titlestr
  oplot,elapsed_time(independent_sample_index)/60,distance(independent_sample_index),psym=1,color=160
  waitkey,dum
end
PRO image_plot,vplot,data_min,data_max,x_range,y_range,ct_string,xlabel,ylabel,title_str


  low_left_corner =[.94,.23]
  sz = [.04,.4]
  img_sz = [.72,.73]
  img_pos = [.1,.15]
  linethick = 1

  ct_start = 0
  ct_stop = 255
  clen = 255
  ct_nix = 0


  data_range = [data_min,data_max]

  image = c_quant(vplot, data_min, data_max, ct_start, ct_stop, ct_nix)

  Setup_vert_ct, low_left_corner, sz, data_range, ct_num, ct_start, ct_stop, ct_string, linethick


!p.title = ''


frame_pos = [img_pos(0), img_pos(1), $
             img_pos(0)+img_sz(0), img_pos(1)+img_sz(1)]

contour, image, /nodata,/noerase, xstyle=1, ystyle=1, $
 color = 1, background = clen-1,xrange=x_range,yrange=y_range,$
 position=frame_pos,/norm,charsize=1.3,xthick=linethick,ythick=linethick,charthick=linethick,xticks=4


imgunder, image

contour, image, /nodata, /noerase, xstyle=1, ystyle=1, $
 xtitle = xlabel, ytitle = ylabel, title = title_str, $
 xrange=x_range,yrange=y_range,position=frame_pos,/norm,$
 color = 255, background = clen-1,charsize=1.3,xthick=linethick,ythick=linethick,charthick=linethick,xticks=4
xmid = mean(x_range)
END

PRO Setup_vert_ct, low_left_corner, sz, data_range, ct_num,ct_start, ct_stop, ct_string,linethick

;print,'low_left_corner, sz, data_range,ct_start, ct_stop, ct_string,linethick',low_left_corner, sz, data_range,ct_start, ct_stop, ct_string,linethick
;clen = ct_len()

pos = [low_left_corner, low_left_corner+sz]

ct = ct_image(ct_start, ct_stop, 'v')

contour, ct, yrange = data_range, $
  /noerase, /nodata, ys=1,xstyle = 5, $
  position = pos, /norm, color = 255, background =0,charsize=1.3,xthick=linethick,ythick=linethick,charthick=linethick
imgunder, ct

contour, ct, yrange = data_range, $
  /noerase, /nodata, ys=1,xstyle = 5, $
  ytitle = ct_string, $
  position = pos, /norm, color = 255, background = 0,charsize=1.3,xthick=linethick,ythick=linethick,charthick=linethick

END

FUNCTION Ct_len

tvlct, r, g, b, /get

sz = size(r)

return, sz(1)
END
FUNCTION C_quant, x, min, max, start, stop, nix

print, 'c_quant: min:', min, ' max:', max, ' start:', start, ' stop:', stop, ' nix:', nix

c = byte(x)

w = where(x LT min, count)

IF (count GT 0) THEN c(w) = byte(nix)

w = where(x GE min, count)

IF (count GT 0) THEN BEGIN

    c(w) = byte((x(w)-min)/(max-min)*(stop-start)+start)

ENDIF

return, c

END


FUNCTION Ct_image, start, stop, dir

IF (dir EQ 'h') THEN BEGIN

    ct = bytarr(stop-start+1, 100)

    FOR i = 0, 99 DO BEGIN

        ct(*, i) = bindgen(stop-start+1)+start

    ENDFOR

ENDIF ELSE IF (dir EQ 'v') THEN BEGIN

    ct = bytarr(100, stop-start+1)

    FOR i = 0, 99 DO BEGIN

        ct(i, *) = bindgen(stop-start+1)+start

    ENDFOR

ENDIF

return, ct

END

    pro imgunder, z, interp=interp, help=hlp

    if (n_params(0) lt 1) or keyword_set(hlp) then begin
      print,' Display image in same area as last plot.'
      print,' imgunder, z'
      print,'   z = scaled byte image to display.  in'
      print,' Keywords:'
      print,'   /INTERP causes bilinear interpolation to be used,'
      print,'     otherwise nearest neighbor interpolation is used.'
      print,' Notes: Do plot,/nodata first to setup plot area,'
      print,'   then use imgunder to display image, then repeat'
      print,'   plot with data, but with /noerase.'
      return
    endif

    if !d.name ne 'PS' then begin
      xx = !x.window*!d.x_size
      yy = !y.window*!d.y_size
      dx = (xx(1) - xx(0))
      dy = (yy(1) - yy(0))
      tv, congrid(z, dx, dy, interp=interp,/minus_one), xx(0), yy(0)+1;add 1 to yy to fix offset of image wrt x axis (12/5/05)
    endif else begin
      xx = !x.window*!d.x_size/!d.x_px_cm
      yy = !y.window*!d.y_size/!d.y_px_cm
      dx = xx(1) - xx(0)
      dy = yy(1) - yy(0)
      tv, z, xx(0), yy(0), xsize=dx, ysize=dy, /cent
    ENDELSE

    return
    END
pro select_index_az,configvars,azimuth,az_proc_index,sweep_count,sweep_count_override,elapsed_time

azmin_proc=configvars.azmin_proc
azmax_proc=configvars.azmax_proc

az_proc_index=where(azimuth gt azmin_proc and azimuth lt azmax_proc,count)

if count eq 0 then begin
   print
   print
   print,'i_az_override = 1'
   print,'processing subset of azimuth angles; no azimuth angles fall in specified range'
   print
   print 
   stop
endif else begin
   wset,0
   plot,elapsed_time/60.,azimuth,ytitle='azimuth angle (deg)',xtitle='time (m)',title='azimuth angles selected for processing (red) ',charsize=1.5,xs=1
   oplot,elapsed_time(az_proc_index)/60.,azimuth(az_proc_index),psym=1,color=160
   sweep_count_override=sweep_count(az_proc_index)
   waitkey,dum
endelse

end;select_index_az
