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

;  Set up arrays to store data

  yy = new ((/nfiles, 500/), "double")
  yy@_FillValue = -999.
  xx = new ((/nfiles, 500/), "double")

;  Intensity corresponding to highest bin to be plotted

  prmax = 200.

  do ifile = 0, (nfiles-1)
    ifile = 0
    f1 = addfile(files(ifile),"r")   ; open file for input

    printVarSummary (f1)

; Note bin_bounds includes both upper and lower limits, so it is dimensioned 1 more than the counts.

    counts_in = f1->num_bin        ;  read in the counts

    printVarSummary (counts_in)

    nbb = dimsizes (counts_in)
    bin_in = f1->bin_bounds(1:nbb)  ;  Read in the bin upper limits. Since bin_bounds should be dimensioned nbb+1, nbb already accounts for the usual n-1

;    Find intensity bin corresponding to upper limit

    npr = closest_val (prmax, bin_in)    
    nb = npr+1

; Select data up to prmax 

    print ("nb, nbb, npr "+nb+" "+nbb+" "+npr)
    y = counts_in(0:nb-1)           ;  select a subset of the input data
    x = bin_in(0:nb-1)

    print (npr)
    printVarSummary (y)
    printVarSummary (x)

;    stop here 27 Mar

;   printVarSummary (tn)
;   printVarSummary (xn)

    printVarSummary (y)
    printVarSummary (yy)

    yy(ifile,0:npr) = y
    xx(ifile,0:npr) = x
    print ("max bin  "+xx(ifile,nb-1))

; xn(j,:) = x(0:nb-1)

; printVarSummary (tn)
; printVarSummary (xn)

;  Convert counts to frequencies

    ktot = f1->numpr
    yy(ifile,:) = yy(ifile,:) / ktot
    yy(ifile,:) = where (yy(ifile,:).gt.0.,yy(ifile,:),yy@_FillValue) 

; do j = 0, (nfiles-1)
;   ktot = sum(tn(j,:))
;   nmeasurable = sum(tn(j,1:nb-1))
;   fmeasurable = int2flt(nmeasurable) / int2flt(ktot)
;   ftrace = int2flt(tn(j,0)) / int2flt(ktot)

;   print ("Number of zero or trace events =    "+tn(j,0))
;   print ("Fraction of zero or trace events = "+ftrace)
;   print ("Total number of measurable precip events = "+nmeasurable)
;   print ("Fraction of measurable precip events     = "+fmeasurable)

; end do

  end do

; print (fn(j,:))



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


  printVarSummary (yy)
  printVarSummary (x)

  xy  = gsn_csm_xy (wks,x,yy,res)                    ; create plot
  draw (xy)

  frame (wks)
 
end

