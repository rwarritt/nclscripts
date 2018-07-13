; ***********************************************
;
; Plot preciptation intensity spectra (pdf)
;
;
; Based on scatter_1.ncl example code
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

  nrows = numAsciiRow ("scatter_in.list")         ;  How many rows in the input file?
  nfiles = nrows - 4                              ;  last four lines are title, marker colors, marker sizes, and marker types
  files = asciiread ("scatter_in.list", nrows, "string")
  print (files)

  rcmname = new ((/100/), "string")

;  Set up arrays to store data, dimensioned (number of files=plots) by (max possible number of bins)

  yy = new ((/nfiles, 500/), "double")
  xx = new ((/nfiles, 500/), "double")
  yy@_FillValue = -999.

;  Intensity corresponding to highest bin to be plotted

  prmax = 200.

;  Loop over the input files

  do ifile = 0, (nfiles-1)
    f1 = addfile(files(ifile),"r")   ; open file for input
    rcmname(ifile) = f1@runname
    printVarSummary (f1)
    counts_in = f1->num_bin          ;  read in the counts
    printVarSummary (counts_in)
    nbb = dimsizes (counts_in)
; Note bin_bounds includes both upper and lower limits, so it is dimensioned 1 more than the counts, i.e. (nbb+1).
    bin_up = f1->bin_bounds(1:nbb)       ;  Read in the bin upper limits. Since bin_bounds should be dimensioned nbb+1, nbb already accounts for the usual n-1
    npr = closest_val (prmax, bin_up)    ;  Find intensity bin corresponding to upper limit
    nb = npr+1                           ;  Account for zero based indexing to get number of bins below upper limit
    print ("nb, nbb, npr "+nb+" "+nbb+" "+npr)
    y = counts_in(0:nb-1)                ; Select data up to prmax 
    x = bin_up(0:nb-1)                   ; x axis is upper limit of this bin
;   printVarSummary (y)
;   printVarSummary (x)
;   printVarSummary (y)
;   printVarSummary (yy)
; Put the x (intensity) and y (count) values into the master arrays
    yy(ifile,0:nb-1) = y
    xx(ifile,0:nb-1) = x
    print ("max bin  "+xx(ifile,nb-1))

; xn(j,:) = x(0:nb-1)
; printVarSummary (tn)
; printVarSummary (xn)

;  Convert counts to frequencies
    ktot = f1->numpr
    yy(ifile,0:nb-1) = yy(ifile,0:nb-1) / ktot

    yy(ifile,0:nb-1) = where (yy(ifile,0:nb-1).gt.0.,yy(ifile,0:nb-1),yy@_FillValue)   ; avoid plotting log of zero

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
; wks   = gsn_open_wks ("ncgm","scatter_rcp85_44")           ; send graphics to NCAR metafile
  wks   = gsn_open_wks ("png","scatter_erai")           ; send graphics to PNG file

;  Set marker types. See https://www.ncl.ucar.edu/Document/Graphics/Images/markers.png
;  Common types:  0 = asterisk,    1 = small filled circle,     2 = plus, 3 = asterisk, 4 = open circle, 5 = cross,
;                 6 = open square, 7 = open triangle, 8 = open inverted triangle, 9 = open diamond, 
;                12 = open five-point star, 16 = large filled cirle

; marker_types = (/16, 6,  6, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16/)

;  Get marker colors, sizes, and types from the last three lines of input data. 
;  Each line is a blank-delimited list.

  mcolors = files(nrows-3)
  msizes = files(nrows-2)
  mtypes = files(nrows-1)

;  Get marker colors, sizes and styles from the last three lines of input.

  delim = " "
  ncolors = str_fields_count (mcolors,delim)
  print ("number of colors = "+ncolors)
  marker_colors = new (ncolors, "string")
  marker_sizes = new (ncolors, "string")
  marker_types = new (ncolors, "string")
  do jc = 0, (ncolors-1)
    marker_colors(jc) = str_get_field (mcolors, (jc+1), delim)
    msize = str_get_field (msizes, (jc+1), delim)
    mtype = str_get_field (mtypes, (jc+1), delim)
    marker_sizes(jc) = stringtofloat (msize)
    marker_types(jc) = stringtointeger (mtype)
    print ("marker_sizes "+marker_sizes(jc)+"  maker_types "+marker_types(jc))
  end do

;  Set up plot resources

  res                   = True                     ; plot mods desired
  res@gsnMaximize       = True                     ; maximize plot

  titlestring = files(nrows-4)
  res@tiMainString      = titlestring              ; add title
  res@tiMainFontHeightF = 0.025                    ; title text size

  res@tiYAxisString    = "Frequency"               ; y-axis label      
  res@tiXAxisString    = "Intensity, mm/d"         ; x-axis label      
  res@tiXAxisFontHeightF = 0.02
  res@tiYAxisFontHeightF = 0.02

  res@xyMarkLineMode    = "Markers"                ; choose to use markers

  res@xyMarkers         =  marker_types(0:nfiles-1)                   ; choose type of marker  
  res@xyMarkerColors := marker_colors(0:nfiles-1)                     ; marker colors
  res@xyMarkerSizes     = marker_sizes(0:nfiles-1)                    ; Marker size (default 0.01)

  res@trYAxisType       = "LogAxis"  ; semilog plot (Y axis log, leave X axis linear)
  res@trXMaxF           = prmax      ;  max value of X axis
  res@trYMaxF           = 1.         ;  max value of Y axis
  res@trYMinF           = 1.e-7      ;  min value of Y axis



  res@vpWidthF      = 0.8                          ;  change the aspect ratio of the plot by scaling the height and width
  res@vpHeightF     = 0.5

  res@gsnFrame               = False               ;  don't advance the frame until we say so

  xy  = gsn_csm_xy (wks,xx(0:nfiles-1,0:nb-1),yy(0:nfiles-1,0:nb-1),res)                    ; create plot

;*****************************************************   
; Manually create and attach legend
;***************************************************** 

  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = 0.0175                 ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification

  res_lines                   = True                  ; polyline mods desired
  res_lines@gsLineDashPattern = 0.                    ; solid line
  res_lines@gsLineThicknessF  = 5.                    ; line thicker
  res_lines@gsLineColor       = "red"                 ; line color

  res_marker = True

  xxl = new ((/nfiles/), "float")
  yyl = new ((/nfiles/), "float")
  dum = new ((/nfiles/), graphic)
  dumt = new ((/nfiles/), graphic)

  yylog = 1.

  do j = 0, (nfiles-1)
    xxl(j) = 0.5 * prmax             ;  x location of legend entry
    yylog = yylog * 0.33             ;  y location of legend entry
    yyl(j) = yylog
    res_marker@gsMarkerIndex = marker_types(j)        
    res_marker@gsMarkerColor = marker_colors(j)
    res_marker@gsMarkerSizeF = marker_sizes(j)
    dum(j) = gsn_add_polymarker(wks,xy,xxl(j),yyl(j),res_marker)              ; add markers
    dumt(j) = gsn_add_text(wks,xy,rcmname(j),(xxl(j)+10.),yyl(j),res_text)    ; add text
  end do


;   All done.  Now draw the plot and advance the frame.

  draw (xy)
  frame (wks)
 
end

