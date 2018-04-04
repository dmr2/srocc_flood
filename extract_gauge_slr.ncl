; get_zos_model_xy.ncl

; written by DJ Rasmussen; last modified Fri May 23 10:37:47 PDT 2014
; 
; Opens gridded sea level file and extracts a specified percentile of 
; sea-level rise for a given lat and lon

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "routines.ncl"
 
 begin

 ; List of tide gauge sites with lat and lon
 lines = asciiread("GPDfits_uhawaii_projectLSLfmt17.tsv",-1,"string")
 nlines = dimsizes(lines)-1   ; First line is a header
 
 delim = "	"; tab delimited
 sid  = str_get_field(lines(1:),2,delim)
 gaugeName = str_get_field(lines(1:),3,delim)
 lat_gauge = tofloat(str_get_field(lines(1:),5,delim))
 lon_gauge = tofloat(str_get_field(lines(1:),6,delim))
 nsite = dimsizes(sid)
 
; convert all site lon to degrees east
 do i=0, nsite-1
  if(lon_gauge(i) .lt. 0 ) then
    lon_gauge(i) = lon_gauge(i) + 360.
  end if
 end do

 filNames = (/"SROCC_RCP26_2050.nc","SROCC_RCP26_2100.nc", \
              "SROCC_RCP45_2050.nc","SROCC_RCP45_2100.nc", \
              "SROCC_RCP85_2050.nc","SROCC_RCP85_2100.nc"/)

 do nf=0,dimsizes(filNames)-1 ; each file

      print("Opening file: "+filNames(nf))
      fil = addfile("./nc/" + filNames(nf),"r")

      lat = fil->y
      lon = fil->x

      slr_mean = fil->total_slr ; (lat,lon)
      slr_std = fil->total_slr_ste ; (lat, lon)

      nlat = dimsizes(dimsizes(lat))

      slr5 = new((/nsite/),"float",-999.999)
      slr50 = new((/nsite/),"float",-999.999)
      slr95 = new((/nsite/),"float",-999.999)

      do i=0, nsite-1

           print("working with site: "+gaugeName(i))
           n = closest_val_AnyOrder(lat_gauge(i),lat)
           m = closest_val_AnyOrder(lon_gauge(i),lon)

           ; If the grid cell nearest the queried lat and lon is land, search the
           ; area for the nearest gridcell that has sea level information
          ; if(sid(i) .eq. 234) then
          ;     print(isnan_ieee(slr_mean(n,m)))
          ;     exit
          ; end if

           if (isnan_ieee(slr_mean(n,m))) then ; try and find nearest grid cell with data

                nval = 5 ; return top 25 nearest grid cells
                ij = getind_n_closest_latlon1d(lat,lon,lat_gauge(i),lon_gauge(i),nval)

                iv = 0
                do while (isnan_ieee(slr_mean(ij(iv,0),ij(iv,1))))

                  if(sid(i) .eq. 234) then
                   ;if(sid(i) .eq. 234) then
                     print(iv+" "+(ij(iv,0)-1)+" "+(ij(iv,1)-1)+" "+slr_mean(ij(iv,0),ij(iv,1)))
                   ;end if
                   print("IN")
                  end if

                   iv=iv+1
                   if(iv.eq.nval*nval) then
                     print("Could not find any SL data nearby.")
                     ;  if(gaugeName(i) .eq. "NEW YORK ( THE BATTERY)") then
                     ;       exit
                     ;  end if
                     iv = iv - 1
                     break
                   end if
                end do

                n = ij(iv,0) ; lat
                m = ij(iv,1) ; lon
                delete(ij)

           end if ; data available in grid cell

      slr5(i) = (slr_mean(n,m) - slr_std(n,m)*1.645 )*100
      slr50(i) = slr_mean(n,m)* 100
      slr95(i) = (slr_mean(n,m) + slr_std(n,m)*1.645 )* 100
      delete([/n,m/])


      end do ; each site

      file_F = str_get_cols(filNames(nf),0,15)+".txt"
      system("/bin/rm -f "+file_F)   ; remove any pre-existing file
      write_table(file_F, "a", [/"Station","ID","Lat","Lon","5th","50th","95th"/], "%s %s %s %s %s %s %s")
      write_table(file_F, "a", [/gaugeName,sid,lat_gauge,lon_gauge,slr5,slr50,slr95/],"%s %s %f %f %f %f %f")


  end do ; each file

end
