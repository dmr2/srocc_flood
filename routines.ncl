
 undef("get_AreaWeight")
 function get_AreaWeight(lat:numeric,lon:numeric)
 local jlat, rad, re, rr, dlon, dx, dy, area
 begin
   jlat = dimsizes( lat )
   rad    = 4.0*atan(1.0)/180.0
   re     = 6371220.0
   rr     = re*rad
 
   dlon   = abs(lon(2)-lon(1))*rr
   dx     = dlon*cos(lat*rad)
   dy     = new ( jlat, typeof(dx))
   dy(0)  = abs(lat(2)-lat(1))*rr
   dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5
   dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr
   area   = dx*dy                                ; cell area function of latitude only
   return(area)
 end

 undef("closest_val_AnyOrder")
 function closest_val_AnyOrder(xVal[1]:numeric, x:numeric)
 local xAbsDif, xMinVal, iClose
 begin
   xAbsDif = abs(xVal-x)
   xMinVal = min(xAbsDif)
   iClose  = (ind(xAbsDif.eq.xMinVal))
   return(iClose(0))   ; original doc says "first occurence"
 end

 undef("getind_n_closest_latlon1d")
 ;
 ; FUNCTION LOCATES NEAREST GRID CELL TO GIVEN LAT AND LON THAT HAS SLR DATA. 
 ;
 function getind_n_closest_latlon1d(lat1d[*]:numeric,lon1d[*]:numeric \
                         ,LAT:numeric, LON:numeric, NVAL:integer)
 local ij, ijdist_nval_lat, ijdist_nval_lon, lat_dist, lon_dist, ilonx, ilatx


 begin

   topNVAL_dist = new((/NVAL*NVAL/),"double",-999.999)
   tmpij = new((/NVAL*NVAL,2/),"integer",-999)
   ij = new((/NVAL*NVAL,2/),"integer",-999)

   r_polar = 6356752.3142
   r_equi = 6378137.0

   lat_dist = 2*(3.141532)*r_polar*(abs(lat1d-LAT)/360.)
   lon_dist = 2*(3.141532)*r_equi*(abs(lon1d-LON)/360.)
  
   idist_ind_lat = toint(dim_pqsort_n(lat_dist,2,0))
   jdist_ind_lon = toint(dim_pqsort_n(lon_dist,2,0))

  ; top NVAL only
   topNVAL_lat = idist_ind_lat(0:NVAL)
   topNVAL_lon = jdist_ind_lon(0:NVAL)

  ; calculate distances using top NVAL 

   icount=0
   do i=0,NVAL-1
    do j=0,NVAL-1

     topNVAL_dist(icount) = sqrt(lat_dist(i)^2 + lon_dist(j)^2)
     tmpij(icount,0) = topNVAL_lat(i)
     tmpij(icount,1) = topNVAL_lon(j)
     icount = icount + 1
    end do
   end do

; return top NVAL*NVAL i,j indices 
   counts = tofloat(dim_pqsort_n(topNVAL_dist,2,0))
   do i=0,NVAL*NVAL-1
     ij(i,0) = tmpij(toint(counts(i)),0)
     ij(i,1) = tmpij(toint(counts(i)),1)
   end do

   ij@long_name = NVAL+" indices that are closest to specified LAT/LON coordinate pairs"

   if (.not.any(ismissing(ij))) then
       delete(ij@_FillValue)
   end if
   return( ij )
 end
