function d=dist_ref(lon_ref,lat_ref,lon_all,lat_all)

dlat=(lat_all-lat_ref).*111.12;
dlon=(lon_all-lon_ref).*(111.12*cosd(lat_ref));

d=sqrt(dlat.^2+dlon.^2);