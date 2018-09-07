;************************************************
begin
;************************************************
;
;  Compute the land-atmosphere coupling metric
;  of Zeng et al. (2010).
;
;  Model output needed:
;   monthly mean precipitation, linked to file name infile_pr.nc
;   monthly mean evapotranspiration, linked to file name infile_ev.nc

;************************************************

;   Start and end years for computations

  year1 = 2076
  yearn = 2095

; Open netCDF files for input:

  f_pr  = addfile("infile_pr.nc","r")     ;  f_pr = monthly precipitation
  f_ev  = addfile("infile_ev.nc","r")     ;  f_ev = monthly evapotranspiration

;   Find out about dimensions.  Ensure that both variables have the same dimensions.

  time_pr = f_pr->time
  time_ev = f_ev->time

  nt_pr = dimsizes (time_pr)
  nt_ev = dimsizes (time_ev)

  if ( dimsizes(time_pr).eq.dimsizes(time_ev) .and. time_pr(0).eq.time_ev(0) .and. time_pr(nt_pr-1).eq.time_ev(nt_ev-1) )  then
    print ("Dimensions match")
    time = time_pr
  else
    print ("Dimensions do not match ")
    printVarSummary (time_pr)
    printVarSummary (time_ev)
  end if

  experiment = f_pr@experiment

  print (experiment)

  lat = f_pr->xlat
  lon = f_pr->xlon

  dimsxy = dimsizes (lat)

  ny = dimsxy(0)
  nx = dimsxy(1)
  nt = dimsizes (time)

  print ("Dimensions nx, ny, nt = "+nx+"  "+ny+"  "+nt)

;   Create a calendar

  cal = cd_calendar (time, -5)


;   Find starting record

  jj = min (ind(cal(:,0).eq.year1)) 

  print ("Starting record for analysis is "+jj+"  "+cal(jj,0)+"  "+cal(jj,1)+"  "+cal(jj,2))


;   Find ending record

  jjn = min (ind(cal(:,0).gt.yearn))
  jjn = jjn - 1
  print ("Ending record for analysis is "+jjn+"  "+cal(jjn,0)+"  "+cal(jjn,1)+"  "+cal(jjn,2))

;   Read the monthly data

  pr = f_pr->pr(jj:jjn,:,:)
  ev = f_ev->evspsbl(jj:jjn,:,:)

;  Convert precipitation and evapotranspiration to mm/day for ease of interpretation.

  pr@_FillValue = default_fillvalue ("float")
  ev@_FillValue = default_fillvalue ("float")

  pr = pr * 86400.
  ev = ev * 86400.

  printVarSummary (pr)
  printVarSummary (ev)

;  Compute monthly climatologies

  pr_anncyc = clmMonTLL (pr)
  ev_anncyc = clmMonTLL (ev)

;  Compute monthly anomalies

  pr_anom = rmMonAnnCycTLL (pr)
  ev_anom = rmMonAnnCycTLL (ev)

;  Compute standard deviations for annual cycle

; pr_std = stdMonTLL (pr_anom)
; ev_std = stdMonTLL (ev_anom)

;  Multiply anomalies to get P'E'

  pprime_eprime = pr_anom * ev_anom

;  Get P' squared

  pprime2 = pr_anom * pr_anom

;  Report number of near-zero values

  pr_min = 0.254                ;  precip rate threshold (mm/day), below which is trace
  ev_min = 0.1                  ;  evapotranspiration threshold
  pp_min = pr_min * pr_min
  pe_min = pr_min * ev_min

  print ("data with pp < "+pp_min+"  mm2/day "+num(pprime2.lt.pp_min))
  print ("data with pp > "+pp_min+"  mm2/day "+num(pprime2.gt.pp_min))
  print ("data with pe < "+pe_min+"   mm2/day "+num(pprime_eprime.lt.pe_min))
  print ("data with pe > "+pe_min+"   mm2/day "+num(pprime_eprime.gt.pe_min))


;  Lazy way to create arrays for the gammas and covariances. 
;  Note "gamma" is an NCL function so we
;  have to use a different name.

  zgamma = pr(0:11,:,:)
  sum_pe = zgamma
  sum_pp = zgamma

  printVarSummary (sum_pp)

;  Compute the gammas

  ndims = dimsizes (pr)
  nt = ndims(0)
  print ("nt "+nt)
  nyears = nt / 12
  print ("nyears "+nyears)
 

;  Initialize the counters

  sum_pe = 0.
  sum_pp = 0.

  print ("Computing covariances")

; Compute the covariances for each month

  do m = 0, 11
    nyr = 0
    j1 = m
    j2 = j1 + (nyears-1) * 12
    do j = j1, j2, 12                 ; sum over years
      sum_pe(m,:,:) = sum_pe(m,:,:) + pprime_eprime(j,:,:)
      sum_pp(m,:,:) = sum_pp(m,:,:) + pprime2(j,:,:)
      nyr = nyr + 1
    end do
    sum_pe(m,:,:) = sum_pe(m,:,:) / nyr
    sum_pp(m,:,:) = sum_pp(m,:,:) / nyr
  end do

  print ("sum_pp < 0.01 "+num(sum_pp.lt.0.01))
  print ("sum_pp > 0.01 "+num(sum_pp.ge.0.01))


;  Avoid analyzing dry areas by putting a floor on the monthly mean precip rate 

; pp_floor = 0.01 * nyears
; pr_min = 0.0001

; sum_pp = where (sum_pp.lt.0.01, sum_pp@_FillValue, sum_pp) 
; printVarSummary (pr_anncyc)
; printVarSummary (sum_pp)

  sum_pp = where (pr_anncyc.lt.pr_min, sum_pp@_FillValue, sum_pp) 

  print ("sum_pp set to missing at "+num(ismissing(sum_pp))+" points")

; print ("sum_pp < 0.01 "+num(sum_pp.lt.0.01))
; print ("sum_pp > 0.01 "+num(sum_pp.ge.0.01))


;  Compute the index

  zgamma = sum_pe / sum_pp

;   write output

  fname = "zgamma.nc"
  f2 = addfile (fname, "c")
  globalAtt = True
  globalAtt@start_year = year1
  globalAtt@end_year   = yearn

  zgamma@long_name = "Zeng gamma"
  sum_pe@long_name = "precipitation evapotranspiration covariance"
  sum_pp@long_name = "precipitation variance"

  zgamma@units = "1"
  sum_pe@units = "1"
  sum_pp@units = "1"

;   These variables have no defined standard name, so delete that attribute

  delete (zgamma@standard_name) 
  delete (sum_pe@standard_name) 
  delete (sum_pp@standard_name) 

  f2->zeng_gamma = zgamma
  f2->pe_covar = sum_pe
  f2->pr_var = sum_pp
  f2->xlat = lat
  f2->xlon = lon

  globalAtt             = True
  globalAtt@start_year = year1
  globalAtt@end_year= yearn
  globalAtt@experiment = experiment
  print (globalAtt@experiment)
  fileattdef (f2, globalAtt) 

end
