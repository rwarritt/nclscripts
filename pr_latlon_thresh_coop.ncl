
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;  
; This script reads daily observed precipitation from Historical
; Climate Network sites and creates frequency distributions of 
; precipitation events falling within a range of continuous "bins."  
; Data are aggregated into a latitude-longitude grid.
;

begin


;
;  *******************************************************************
;
;                     USER SETTABLE PARAMETERS
;
;  *******************************************************************
;

; Max and min lat/lon

latlon = asciiread ("latlon.in", -1, "float")

latmin = latlon(0)
latmax = latlon(1)
lonmin = latlon(2)
lonmax = latlon(3)

; Ensure lons are degreese east

print ("latmin, latmax, lonmin, lonmax "+latmin+" "+latmax+" "+lonmin+" "+lonmax)
lonmin = mod ((lonmin+360.), 360.)
lonmax = mod ((lonmax+360.), 360.)

; check if lat/lon are flipped

if (lonmin .gt. lonmax) then
  ltmp = lonmin
  lonmin = lonmax
  lonmax = ltmp
end if

print ("latmin, latmax, lonmin, lonmax "+latmin+" "+latmax+" "+lonmin+" "+lonmax)

; Start and end year, and season

; for ERAI runs
;start_year = 1990
;end_year = 2009

; for GCM runs
start_year = 1981
end_year = 2005

select_season = "warm"


; Get the list of stations

nfiles = numAsciiRow ("coop/allstations.txt")         ;  How many files?

print ("Preparing to read list of "+nfiles+" files")

files = asciiread ("coop/allstations.txt", nfiles, "string")

printVarSummary (files)

; Set minimum value for measurable precip (mm/day). Below this is considered trace.

trace = 0.1                    
trace@long_name = "precip_lower_limit"
trace@units = "mm d-1"

; Set maximum fraction of missing data for a station to be included

maxmiss = 0.1              

do_exceed = False
do_pdf_region = True


; *****************   
; *****************     end user settable parameters
; *****************  

; Check that the season is an allowable value

if (select_season.ne."warm" .and.  select_season.ne."cool" .and.  select_season.ne."year") then
  print ("STOP - Unknown season "+select_season)
end if

; Open file for output

f2 = addfile ("binned_coop.nc", "c")

; Set up bins for pdf

  opt = True
  opt@bin_spacing = 2.
  opt@bin_min =  0.
  opt@bin_max =  500.

  xbins = (opt@bin_max - opt@bin_min) / opt@bin_spacing
  print ("xbins = "+xbins)

  nbins = round (xbins, 3)

  num_bin = new (nbins, double)
  num_bin!0 = "nbins"
  num_bin@long_name = "number of events in each intensity bin"
  num_bin = 0

  bin_center = new (nbins, double)
  bin_center!0 = "nbins"
  bin_center@long_name = "center value of intensity bin"
  bin_center@units = "mm d-1"
  bin_center = 0

  bin_bounds = new (nbins+1, double)
  bin_bounds!0 = "nbounds"
  bin_bounds@long_name = "bounds of intensity bin"
  bin_bounds@comment = "lowest bound should equal trace"
  bin_bounds@units = "mm d-1"
  bin_bounds = 0

  printVarSummary (num_bin)
  printVarSummary (bin_center)

  nmax_tot = 0
  nuse_tot = 0
  nlo_tot = 0
  nhi_tot = 0
  numt_tot = 0
  numpr_tot = 0
  numtrace_tot = 0
  numpoints = 0

; Find first and last records to be used

  fin = files(0)
  f0 = addfile (fin, "r")
  timein = f0->time
  nt = dimsizes (timein)

; Turn the time into an understandable calendar date

  cal = cd_calendar (timein, -5)  ; note year = calendar(:,0) and month is calendar (:,1)

  tind = ind (cal(:,0).ge.start_year .and. cal(:,0).le.end_year .and. cal(:,1).ge.4 .and. cal(:,1).le.9) 
  nti = dimsizes (tind)

  tcal = cal(tind,:)

; do j = 0, (nti-1)
;   print (tcal(j,0)+" "+tcal(j,1)+" "+tcal(j,2))
; end do


; itstart is the first record in the starting year.
; Report the year corresponding to this record.

  istarty = get1Dindex (cal(:,0), start_year)
  print ("Starting record: istarty, year = "+istarty+"   "+cal(istarty,0))

  print ("Season "+select_season)
  print ("istarty "+istarty)

  if (select_season .eq. "warm") then
    istartm = get1Dindex (cal(istarty:,1), 4)
    istartm = istartm + istarty
  end if

  print ("Starting record: istartm, month = "+istartm+"   "+cal(istartm,0))

;  trust but verify

  if (istarty .gt. 0) then
    print ("Year for record (istarty-1) = "+cal(istarty-1,0))
  end if

  if (istartm .gt. 0) then
    print ("Year for record (istartm) = "+cal(istartm,0))
    print ("Year for record (istartm-1) = "+cal(istartm-1,0))
    print ("Month for record (istartm)  = "+cal(istartm,1))
    print ("Month for record (istartm-1)  = "+cal(istartm-1,1))
  end if

; Find the ending record

  iendy = nt - 1

  if (cal(nt-1,0) .gt. end_year) then
    iendy = max (ind(cal(:,0).eq.end_year))           ;  index of last time within the ending year
  end if

  print ("Ending record: iendy, year = "+iendy+"   "+cal(iendy,0)+"  "+cal(iendy,1)+"  "+cal(iendy,2))

  if (iendy .lt. (nt-1)) then
    print ("Year for record (iendy+1) = "+cal(iendy+1,0))
  end if


  if (select_season .eq. "warm") then
    iendm = max (ind ((cal(:,0).eq.end_year).and.(cal(:,1).eq.9)))            
  end if

  print ("iendm "+iendm+" "+cal(iendm,0)+" "+cal(iendm,1)+" "+cal(iendm,2))

  if (iendm .lt. (nt-1)) then
    print ("Year for record (iendm) = "+cal(iendm,0))
    print ("Year for record (iendm-1) = "+cal(iendm-1,0))
    print ("Month for record (iendm)  = "+cal(iendm,1))
    print ("Month for record (iendm+1)  = "+cal(iendm+1,1))
  end if

  itstart = istartm
  itend = iendm
  print ("Starting record: itstart, month, year = "+itstart+"  "+cal(itstart,1)+"  "+cal(itstart,0)+"  "+cal(itstart,2))
  print ("Ending record: itend, month, year     = "+itend+"  "+cal(itend,1)+"  "+cal(itend,0)+"  "+cal(itend,2))

;  Loop over stations

if (do_pdf_region) then

  points_lat = new (50000, typeof(f0->lat))
  points_lon = new (50000, typeof(f0->lon))

;  npdf = number of events in each bin for this grid point
;  numpdf = total number of events summed over all bins for this grid point
;  num_bin = number of events in each bin accumulated over all grid points
;  num_bin_tot = total number of events summed over all bins and grid points
;  nmax_tot = total number of non-missing (i.e., non-trace) values, accumulated over all grid points
;  nuse_tot = total number of values used (i.e., lying within the bin limits), accumulated over all grid points
;  nlo_tot = total number of values below the lowest bin, accumulated over all grid points (should always be zero)
;  nhi_tot = total number of values above the highest bin, accumulated over all grid points 
;  npoints = number of locations (stations) used


  do i = 0, (nfiles-1)
    flag = 0
    fin = files(i)
    f1 = addfile (fin, "r")
    lat = f1->lat
    lon = f1->lon
    lon = mod (lon + 360., 360.)          ;    convert to degrees east
; check if lat is within our region
    if ((lat.lt.latmin) .or. (lat.gt.latmax)) then
      print (fin+"  not in target region")
      continue
    end if
; check if lon is within our region
    if ((lon.lt.lonmin) .or. (lon.gt.lonmax)) then
      print (fin+"  not in target region "+lat+"  "+lon)
      continue
    end if
    print (fin+"  is in target region")
    time = f1->time
    nt = dimsizes (time)
    cale = cd_calendar (time, -5)
    pr_in = f1->pr 
; use only warm season
    if (select_season .eq. "warm") then
      prind = ind (cale(:,0).ge.start_year .and. cale(:,0).le.end_year .and. cale(:,1).ge.4 .and. cale(:,1).le.9) 
      pr = pr_in(prind)
    end if
; check data continuity
    nmiss = num(ismissing(pr))
    ntot = dimsizes (pr)
    print ("nmiss, ntot "+nmiss+" "+ntot)
    missfrac = nmiss / int2flt(ntot)
    print ("fraction of missing values = "+missfrac)
; bail out if fraction of missing data is too large
    if (missfrac .gt. maxmiss) then
      print ("too many missing values, nmiss, ntot, nfrac = "+nmiss+" "+ntot+" "+missfrac)
      continue
    end if
;
; at this point we have a station within our target region that has the required amount of good data
;
    flag = 1 
    numpoints = numpoints + 1                    ;  increment number of locations in our sample
    points_lat(numpoints-1) = (/ lat /)
    points_lon(numpoints-1) = (/ lon /)
    pr = pr * 86400.                             ;  convert precip to mm/day
    pr = where (pr.gt.trace, pr, pr@_FillValue)  ;  ignore values below cutoff (trace) by setting them to missing
    numpr = num (.not.ismissing(pr))             ;  number of non-trace values at this grid point
    numtrace = num (ismissing(pr))               ;  number of trace values at this grid point
    numtrace = numtrace - nmiss                  ;  adjust for the number of obs that really are missing
    numt = dimsizes(pr)                          ;  number of values at this grid point, both trace and non-trace (equals number of time steps)
    numpr_tot = numpr_tot + numpr                ;  accumulate number of non-trace values over all grid points
    numtrace_tot = numtrace_tot + numtrace       ;  accumulate number of trace values
    numt_tot = numt_tot + numt                   ;  accumulate total number of values over all grid points
    pr_pdf = pdfx(pr, nbins, opt)                 ;  compute pdf here -- note pdfx computes pdf as percent
    npdf = pr_pdf@nMax * pr_pdf / 100            ;  convert percent to counts
    nmax_tot = nmax_tot + pr_pdf@nMax            ;  nMax is how many non-missing values were input
    nuse_tot = nuse_tot + pr_pdf@nUse            ;  nUse is how many values were used in the pdf, i.e., within the bin ranges
    nehi = 0
    nelo = 0
    if (isatt(pr_pdf,"nLoOut")) then
      nelo = pr_pdf@nLoOut                 ;  number of values at this grid point that are below the lowest bin
    end if
    nlo_tot = nlo_tot + nelo               ;  accumulate number of values below lowest bin
    if (isatt(pr_pdf,"nHiOut")) then
      nehi = pr_pdf@nHiOut
    end if
    nhi_tot = nhi_tot + nehi               ;  number of values at this grid point that are above highest bin
    num_bin = num_bin + npdf               ;  add counts to accumulator bins
    pr_max = max (pr)
    num_bin_tot = sum (num_bin)
    numpdf = sum (npdf)
    print ("numpdf, num_bin_tot, nmax_tot, nuse_tot, nlo_tot, nhi_tot "+numpdf+"  "+num_bin_tot+"  "+nmax_tot+"  "+nuse_tot+"  "+nlo_tot+"  "+nhi_tot+"  "+pr_max+" "+nehi)
;   print (pr_pdf)
;   print (npdf)
;   print (num_bin)

    if (numpoints .eq. 1) then
      bin_center = pr_pdf@bin_center
      bin_bounds = pr_pdf@bin_bounds
    end if
    if (flag .eq. 1) then
      delete (pr)
      delete (pr_pdf)       ;     for good measure
    end if
    delete (f1)

;   if (numpoints .gt. 3) then
;     break                  ;  bail out of the loop
;   end if

;   exit

  end do       ; end loop over stations

end if         ;  end if for do_pdf_region

  print ("numpoints  "+numpoints)

  printVarSummary (numpoints)

  do jpts = 0, (numpoints-1)
    print ("lat, lon  "+points_lat(jpts)+"  "+points_lon(jpts))
  end do

  print ("num_bin_tot, nmax_tot, nuse_tot, nlo_tot, nhi_tot "+num_bin_tot+"  "+nmax_tot+"  "+nuse_tot+"  "+nlo_tot+"  "+nhi_tot)
  print ("numt_tot, numpr_tot, numtrace_tot  "+numt_tot+"  "+numpr_tot+"  "+numtrace_tot)

  printVarSummary (points_lat)
  printVarSummary (points_lon)

  lat_points = new (numpoints, typeof(points_lat))
  lon_points = new (numpoints, typeof(points_lon))

  lat_points = points_lat(0:numpoints-1)
  lon_points = points_lon(0:numpoints-1)
  lat_points!0 = "numpoints"
  lon_points!0 = "numpoints"
  lat_points@long_name = "latitude of data points used"
  lon_points@long_name = "longitude of data points used"
  lat_points@units = "degrees north"
  lon_points@units = "degrees east"

  bin_center!0 = "nbins"

  printVarSummary (numpoints)

  f2->num_bin = num_bin
  f2->lat_points = lat_points
  f2->lon_points = lon_points
  f2->bin_center = bin_center
  f2->bin_bounds = bin_bounds

  print (trace)
  trace!0 = "ncl_scalar"
  print (trace)

  numt_tot!0 = "ncl_scalar"
  numpr_tot!0 = "ncl_scalar"
  numtrace_tot!0 = "ncl_scalar"

  latmin!0 = "ncl_scalar"
  latmax!0 = "ncl_scalar"
  lonmin!0 = "ncl_scalar"
  lonmax!0 = "ncl_scalar"

  latmin@long_name = "south limit of latitude"
  latmax@long_name = "north limit of latitude"
  lonmin@long_name = "west limit of longitude"
  lonmax@long_name = "east limit of longitude"

  latmin@units = "degrees north"
  latmax@units = "degrees north"
  lonmin@units = "degrees east"
  lonmax@units = "degrees east"

  numt_tot@long_name = "total number of data values"
  numpr_tot@long_name = "number of values with precipitation trace or above"
  numtrace_tot@long_name = "number of values with precipitation less than trace"

  f2->trace = trace
  f2->numt = numt_tot                  
  f2->numpr = numpr_tot                 ;  accumulate number of non-trace values
  f2->numtrace = numtrace_tot
  f2->latmin = latmin
  f2->latmax = latmax
  f2->lonmin = lonmin
  f2->lonmax = lonmax

; copy_VarAtts(f1, f2)             ;   copy global attributes from input to output

;  more global attributes
  globalAtt = True
  globalAtt@start_year = start_year
  globalAtt@end_year = end_year
  fileattdef (f2, globalAtt)
 
  exit 

;         num_bins(:,jj,ii) = num_bins(:,jj,ii) + npdf
;         print ("numx, numt = "+numx+"  "+numt)

;------------------------------------------------------------


  print ("Done.")


  printVarSummary (lat)
  printVarSummary (lon)

; f2->lat = lat
; f2->lon = lon



end
