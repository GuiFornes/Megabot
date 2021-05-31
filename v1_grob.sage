R.<dx,dy,hx,hy,ix,iy,ex,ey,fx,fy,gx,gy,jx,jy,ao,bc,bf,ae,de,ef,fg,fh,gi,gj> = QQ[]
equations=[]
equations.append(ao-0.150)
equations.append(bc-0.300)
equations.append(bf-0.600)
equations.append(ae-0.500)
equations.append(de-0.100)
equations.append(ef-0.450)
equations.append(fg-0.300)
equations.append(fh-0.300)
equations.append(gi-0.500)
equations.append(gj-1.000)
equations.append(jx-0.79)
equations.append(jy+0.49)

equations.append( (ex-0.15)^2 + ey^2 - ae^2 )
equations.append( (fx-ex)^2 + (fy-ey)^2 - ef^2 )
equations.append( (fx-gx)^2 + (fy-gy)^2 - fg^2 )
equations.append( (gx-ex)^2 + (gy-ey)^2 - (fg+ef)^2 )

equations.append( (fx-hx)^2 + (fy-hy)^2 - fh^2 )
equations.append( (fx-0.1)^2 + (fy)^2 - bf^2 )
equations.append( (hx-0.1)^2 + (hy)^2 - (bf-fh)^2 )

equations.append( (dx-0.15)^2 + dy^2 - (ae-de)^2)
equations.append( (dx-ex)^2 + (dy-ey)^2 - (de)^2)

equations.append( (gx-ix)^2 + (gy-iy)^2 - (gi)^2)
equations.append( (gx-jx)^2 + (gy-jy)^2 - (gj)^2)
equations.append( (jx-ix)^2 + (jy-iy)^2 - (gj+gi)^2)
I = R*list(equations)
I.groebner_basis()
I.dimension()