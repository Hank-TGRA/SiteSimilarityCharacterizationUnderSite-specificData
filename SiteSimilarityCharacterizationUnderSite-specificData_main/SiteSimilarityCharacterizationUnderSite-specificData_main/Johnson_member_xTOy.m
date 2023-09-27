function yi = Johnson_member_xTOy(x, jmtype, jparas)

switch jmtype
    case 'SB'
        ax = jparas(1); bx = jparas(2); ay = jparas(3); by = jparas(4);
        xni = (x-bx)./ax;
        intervar = exp(xni);
        yni = (intervar./(1+intervar));
        yi = by+ay*yni;
    case 'SU'
        ax = jparas(1); bx = jparas(2); ay = jparas(3); by = jparas(4);
        xni = (x-bx)./ax;
        yni = sinh(xni);
        yi = by+ay*yni;
    case 'SL'
        ax = jparas(1); bx2 = jparas(2); by = jparas(3);
        xni = (x-bx2)./ax;
        yni = exp(xni);
        yi = by+yni;
end
end