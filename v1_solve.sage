bx,by,jx,jy,fx,fy,gx,gy,bf,fg,gj=var('bx by jx jy fx fy gx gy bf fg gj')

eq_dist_bf=(fx-bx)^2 + (fy-by)^2 == bf^2
eq_dist_fg=(fx-gx)^2 + (fy-gy)^2 == fg^2
eq_dist_gj=(gx-jx)^2 + (gy-jy)^2 == gj^2
eq_pytha = bf^2 + fg^2 + gj^2 + (jx -bx)^2 + (jy - by)^2 == (fx -jx)^2 + (fy-jy)^2 + (gx -bx)^2 + (gy-by)^2 + 4 * ( ((fx+jx)/2 - (gx+bx)/2)^2 + ((fy+jy)/2 - (gy+by)/2)^2 )
solve([eq_dist_bf,eq_dist_fg,eq_dist_gj,eq_pytha],fx,fy,gx,gy)



jx,fx,fy,gx,gy,bf,fg,gj=var('jx fx fy gx gy bf fg gj')

eq_dist_bf=(fx)^2 + (fy)^2 == bf^2
eq_dist_fg=(fx-gx)^2 + (fy-gy)^2 == fg^2
eq_dist_gj=(gx-jx)^2 + (gy)^2 == gj^2
eq_pytha = bf^2 + fg^2 + gj^2 + (jx)^2 == (fx -jx)^2 + (fy)^2 + (gx)^2 + (gy)^2 + 4 * ( ((fx+jx)/2 - (gx)/2)^2 + ((fy)/2 - (gy)/2)^2 )
solve([eq_dist_bf,eq_dist_fg,eq_dist_gj,eq_pytha],fx,fy,gx,gy)


