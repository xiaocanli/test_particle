;******************************************************************************
; Procedures to deal with velocity fields.
;******************************************************************************
PRO velocity
    common dims, nx, ny, nz, nt
    nx = 256
    ny = 256
    nz = 256
    nt = 16
    print, "Procedures to deal with velocity fields."
    dz = 16
    nz = 256
    it = 8
    for iz = 0, nz, dz do begin
        PlotV2d, it, iz
    endfor
END

;******************************************************************************
; Read 4-dimensional velocity field data from HDF5 file.
; Input:
;   fname, gname, dname: file name, group name and data set name.
;   count, offset: the dimensions and offset in the whole data set for
;       hyperslab selection.
; Return:
;   data: the data set readed from the file.
;******************************************************************************
pro ReadDataHDF5, fname, gname, dname, count, offset, data
    ; Open the HDF5 file.
    file_id = H5F_OPEN(fname)
    ; Open the group and dataset
    group_id = H5G_OPEN(file_id, gname)
    dset_id = H5D_OPEN(group_id, dname)
    ; Open up the dataspace associated wit the dataset
    filespace = H5D_GET_SPACE(dset_id)
    ;dims = H5S_GET_SIMPLE_EXTENT_DIMS(filespace)
    H5S_SELECT_HYPERSLAB, filespace, offset, count, /RESET
    ; Create a simple dataspace to hold the result.
    memspace = H5S_CREATE_SIMPLE(count)
    ; Read in the actual data.
    data = H5D_READ(dset_id, FILE_SPACE=filespace, $
        MEMORY_SPACE=memspace)
    H5S_CLOSE, memspace
    H5S_CLOSE, filespace
    H5D_CLOSE, dset_id
    H5G_CLOSE, group_id
    H5F_CLOSE, file_id
end

;******************************************************************************
; Plot 2D contour of velocity fields
; Input:
;   it: the ID of time slice.
;   iz: the ID of the z-direction slice.
;******************************************************************************
pro PlotV2d, it, iz
    common dims
    offset = [0, 0, iz, it]
    count = [nx, ny, 1, 1]
    ReadDataHDF5, 'u4d.h5', '/u4d', 'uz4d', count, offset, data

    x = findgen(nx)
    y = findgen(ny)
    data = reform(data)
    data1 = smooth(data, [1,1])
    maxd = max(data1)
    mind = min(data1)
    lims1 = [mind, maxd]
    BlueWhiteRed, rgbtable, lims1
    title = 'z = ' + string(iz, format='(I3.3)')
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
    fname = 'v2d_xy_' + string(iz, format='(I3.3)') + '.jpg'
    im1.save, fname
    im1.close
end

;******************************************************************************
; Plot the time evolution of velocity field at one point.
; Input:
;   ix, iy, iz: the position of the chosen point.
;******************************************************************************
pro PlotVelocityEvolution, ix, iy, iz
    common dims
    offset = [ix, iy, iz, 0]
    count = [1, 1, 1, nt]
    dim_id = 'z'
    dname = 'u' + dim_id + '4d'
    ReadDataHDF5, 'u4d.h5', '/u4d', dname, count, offset, data
    t = findgen(nt)
    ytitle = '$u_' + dim_id + '$'
    p1 = plot(t, data, 'k2', $
        xtitle='t', ytitle=ytitle, $
        font_size=24)
end

;******************************************************************************
; Description: Blue-White-Red color table.
;   Originally from matlab file exchange: 
;   http://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
;******************************************************************************
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

; Routine to create a video
PRO Video_current, qname, nt, icumulation
    width = 768
    height = 360
    IF (icumulation EQ 1) THEN BEGIN
        height = 480
    ENDIF
    frames = nt
    fps = 10

    ; Create object and initialize video/audio streams
    fname = qname + '.mp4'
    oVid = IDLffVideoWrite(fname)
    vidStream = oVid.AddVideoStream(width, height, fps)
   
    ; Generate video frames
    FOR i = 0, frames-1 DO BEGIN
        print, i
;        lims = [-1.0,1.0]
;        firehose, i, 70, lims, im
        lims = [-0.0005,0.0005]
        testplot, qname, i, icolor, lims, im
;        lims = 5.0
;        preanisotropy, i, 70, lims, im
;        fname = 'img/jy_all00' + STRING(i+1,FORMAT='(I2.2)') + '.png'
;        im = IMAGE(/buffer,fname,dimensions=[width,height])
        time = oVid.Put(vidStream, im.CopyWindow())
        im.Close
    ENDFOR
    ;im.Close
    
    ; Close the file
    oVid.Cleanup
END
