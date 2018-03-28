; ***********************************************
;
; Scatter plot with semilog axes and multiple data sets.
;
;
; Based on scatter_1.ncl example code
;
; Concepts illustrated:
;   - Drawing a scatter plot
;   - Changing the markers in an XY plot
;   - Changing the marker color in an XY plot
;   - Changing the marker size in an XY plot
;
; ***********************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin

;************************************************
; Open data file
;************************************************


;  Read the list of file names

  nfiles = numAsciiRow ("scatter_in.list")         ;  How many files?
  files = asciiread ("scatter_in.list", nfiles, "string")

  print (files)

  ifile = 0
  f1 = addfile(files(ifile),"r")

  printVarSummary (f1)

  tin = f1->num_bin
  nbb = dimsizes (tin)
  xin = f1->bin_bounds
  xx = xin(1:nbb)

;  Find intensity bin corresponding to upper limit

  prmax = 200.
  npr = closest_val (prmax, xx)
  nb = npr+1

  t = tin(0:npr)
  x = xx(0:npr)

  printVarSummary (t)
  printVarSummary (x)


  nfiles = 1


  tn = new ((/nfiles, nb/), "double")
  xn = new ((/nfiles, nb/), "double")
  fn = tn

  printVarSummary (tn)
  printVarSummary (xn)


  j = 0

  tn(j,:) = t
  xn(j,:) = x(0:nb-1)

; xn(j,:) = x(0:nb-1)

  printVarSummary (tn)
  printVarSummary (xn)

  print ("max bin  "+xn(j,nb-1))

  data2 = new ((/2,nb/), "double")
  data2@_FillValue = -999.

  ktot = f1->numpr


;  Convert counts to frequencies

; do j = 0, (nfiles-1)
;   ktot = sum(tn(j,:))
;   nmeasurable = sum(tn(j,1:nb-1))
;   fmeasurable = int2flt(nmeasurable) / int2flt(ktot)
;   ftrace = int2flt(tn(j,0)) / int2flt(ktot)

;   print ("Number of zero or trace events =    "+tn(j,0))
;   print ("Fraction of zero or trace events = "+ftrace)
;   print ("Total number of measurable precip events = "+nmeasurable)
;   print ("Fraction of measurable precip events     = "+fmeasurable)

    fn(j,:) = tn(j,:)
    fn(j,:) = fn(j,:) / ktot
    data2(j,:) = fn(j,:)
    data2(j,:) = where (data2(j,:).gt.0.,data2(j,:),data2@_FillValue) 
; end do



  print (fn(j,:))



;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("ncgm","scatter")           ; send graphics to PNG file

;  Set marker types. See https://www.ncl.ucar.edu/Document/Graphics/Images/markers.png
;  Common types:  0 = asterisk,    1 = small filled circle,     2 = plus, 3 = asterisk, 4 = open circle, 5 = cross,
;                 6 = open square, 7 = open triangle, 8 = open inverted triangle, 9 = open diamond, 
;                12 = open five-point star, 16 = large filled cirle

  marker_types = (/16, 6, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16/)
  marker_colors = (/"blue", "black", "green", "black", "blue", "red", "green", "black", "blue", "red", "green", "black", "blue", "red", "green", "black" /)
  marker_sizes = (/ 0.002, 0.004, 0.002, 0.004, 0.002, 0.002, 0.004, 0.004, 0.002, 0.002, 0.004, 0.004, 0.002, 0.002, 0.004, 0.004 /)

  res                   = True                     ; plot mods desired
  res@gsnMaximize       = True                     ; maximize plot
  res@tiMainString      = "Scatter Plot"           ; add title
  res@xyMarkLineMode    = "Markers"                ; choose to use markers

  res@xyMarkers         =  marker_types(0:nfiles-1)                     ; choose type of marker  
  res@xyMarkerColors := marker_colors(0:nfiles-1)   ; 2 different colors
  res@xyMarkerSizes     = marker_sizes(0:nfiles-1)                    ; Marker size (default 0.01)

  res@trYAxisType       = "LogAxis"  ; semilog plot (Y axis log, leave X axis linear)
  res@trXMaxF           = 200.       ;  max value of X axis
  res@trYMaxF           = 1.         ;  max value of Y axis
  res@trYMinF           = 1.e-7      ;  min value of Y axis

;   Use these variables to change the aspect ratio of the plot by scaling the height and width

  res@vpWidthF      = 0.8
  res@vpHeightF     = 0.5

; res@gsnYAxisIrregular2Log = True
; plot  = gsn_csm_y (wks,t,res)                    ; create plot
  res@gsnFrame               = False


  printVarSummary (data2)
  printVarSummary (x)

  xy  = gsn_csm_xy (wks,x,data2,res)                    ; create plot
  draw (xy)

  frame (wks)
 
end

