pro bfield
    print, "Procedures to plot magnetic field."
    common dims, nmax, dxyz, lmax
    nmax = 256
    lmax = 20.0
    dxyz = lmax / nmax
    bx = dblarr(nmax, nmax)
    by = dblarr(nmax, nmax)
    bz = dblarr(nmax, nmax)
    for iz = 0, nmax-1, 16 do begin
        GetBfieldArray, iz, bx, by, bz
        PlotB2d, iz, bx, 'x'
        PlotB2d, iz, by, 'y'
        PlotB2d, iz, bz, 'z'
    endfor
end

pro GetBfieldArray, iz, bx, by, bz
    common dims
    z = dxyz * iz
    for i = 0, nmax-1 do begin
        x = i * dxyz
        for j = 0, nmax-1 do begin
            y = j * dxyz
            ForceFreeField, x, y, z, $
                bx1, by1, bz1
            bx[i,j] = bx1
            by[i,j] = by1
            bz[i,j] = bz1
        endfor
    endfor
end

;******************************************************************************
; Procedure to calculate force-free field
; Input:
;   x, y, z: spatial positions.
; Return:
;   bx, by, bz: the resulting magnetic field.
;******************************************************************************
pro ForceFreeField, x, y, z, bx, by, bz
    ffA = 1.0
    ffB = sqrt(2.0/3.0)
    ffC = sqrt(1.0/2.0)
    bx = ffA*sin(z) + ffC*cos(y);
    by = ffB*sin(x) + ffA*cos(z);
    bz = ffC*sin(y) + ffB*cos(x);
end

pro PlotB2d, iz, data, qname
    common dims
    x = findgen(nmax)
    y = findgen(nmax)
    data = reform(data)
    data1 = smooth(data, [1,1])
    maxd = max(data1)
    mind = min(data1)
    lims1 = [mind, maxd]
    BlueWhiteRed, rgbtable, lims1
    title = '$B_' + qname + '$' + ', z = ' $
        + string(iz, format='(I3.3)')
    im1 = image(data1, x, y, /buffer, $
        xtitle = 'x', ytitle = 'y', $
        title = title, $
        font_size=16, $
        position=[0.10, 0.15, 0.85, 0.9], $
        rgb_table=rgbtable, axis_style=2, interpolate=1)

    print, 'Maximum', maxd
    print, 'Minimum', mind
    im1.max_value = lims1(1)
    im1.min_value = lims1(0)
    ;pos1 = im1.position
    ;print, pos1
    CB1 = COLORBAR(TARGET=im1,ORIENTATION=1,$
        position=[0.83,0.15,0.86,0.9])
    CB1.TEXTPOS = 1
    CB1.TICKDIR = 1
    CB1.FONT_SIZE = 16
    ;CB1.TITLE='$|B|$'
    fname = 'b' + qname + '2d_xy_' + string(iz, format='(I3.3)') + '.jpg'
    im1.save, fname
    im1.close
end

PRO BlueWhiteRed, rgbtable, lims
    p = FLTARR(3,5)
    p = [[0.0,0.0,0.5],[0.0,0.5,1.0],[1.0,1.0,1.0],[1.0,0.0,0.0],[0.5,0.0,0.0]]
    ;p = [[0.0,0.0,0.5],[0.0,0.5,1.0],[0.5,0.5,0.5],[1.0,0.0,0.0],[0.5,0.0,0.0]]
    m = 256
    rgbtable = FLTARR(3,m)
    IF ((lims(0) LT 0) AND (lims(1) GT 0)) THEN BEGIN
      ; It has both negative and positive values.
      ; Find portion of negative values.
        ratio = abs(lims(0)) / (abs(lims(0))+abs(lims(1)))
        neglen = round(m*ratio)
        poslen = m - neglen

      ; Colorbar for negative values
        cindex = FINDGEN(neglen)*2.0/neglen
        rgbtable(0,0:neglen-1) = INTERPOLATE(p(0,0:2),cindex)
        rgbtable(1,0:neglen-1) = INTERPOLATE(p(1,0:2),cindex)
        rgbtable(2,0:neglen-1) = INTERPOLATE(p(2,0:2),cindex)
        rgbtable(*,0:neglen-1) = rgbtable(*,0:neglen-1)*255

      ; Colorbar for positive values
        cindex = FINDGEN(poslen)*2.0/poslen
        rgbtable(0,neglen:m-1) = INTERPOLATE(p(0,2:4),cindex)
        rgbtable(1,neglen:m-1) = INTERPOLATE(p(1,2:4),cindex)
        rgbtable(2,neglen:m-1) = INTERPOLATE(p(2,2:4),cindex)
        rgbtable(*,neglen:m-1) = rgbtable(*,neglen:m-1)*255
    ENDIF ELSE IF (lims(0) GT 0) THEN BEGIN
      ; Just positive values
        cindex = FINDGEN(m)*2.0/m
        rgbtable(0,*) = INTERPOLATE(p(0,2:4),cindex)
        rgbtable(1,*) = INTERPOLATE(p(1,2:4),cindex)
        rgbtable(2,*) = INTERPOLATE(p(2,2:4),cindex)
        rgbtable = rgbtable*255
    ENDIF ELSE IF (lims(1) LT 0) THEN BEGIN
      ; Just negative vaules
        cindex = FINDGEN(m)*2.0/m
        rgbtable(0,*) = INTERPOLATE(p(0,2:0:-1),cindex)
        rgbtable(1,*) = INTERPOLATE(p(1,2:0:-1),cindex)
        rgbtable(2,*) = INTERPOLATE(p(2,2:0:-1),cindex)
        rgbtable = rgbtable*255
    ENDIF
END

pro PlotBdistribution
    fname = 'bdist.dat'
    nlines = FILE_LINES(fname)
    data = dblarr(2,nlines)
    openr, lun, fname, /get_lun
    readf, lun, data
    free_lun, lun 

    p1 = plot(data(0,*), data(1,*), $
        xtitle='$|B|$', ytitle='$dN/d|B|$', $ 
        thick=2, /ylog, font_size=24, $
        position=[0.2,0.2,0.95,0.95])
    p1.save, 'img/b_dist.eps'
end
