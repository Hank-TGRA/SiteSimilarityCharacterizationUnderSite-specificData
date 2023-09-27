function x = Johnson_member_yTOx(y, jmtype, jparas)
switch jmtype
    case 'SB'
        ax = jparas(1); bx = jparas(2); ay = jparas(3); by = jparas(4);
        yn = (y-by)./ay;
        x = bx+ax.*log(yn./(1-yn));
    case 'SU'
        ax = jparas(1); bx = jparas(2); ay = jparas(3); by = jparas(4);
        yn = (y-by)./ay;
        x = bx+ax.*log(yn+sqrt(1+yn.^2));
    case 'SL'
        ax = jparas(1); bx2 = jparas(2); by = jparas(3);
        x = bx2+ax.*log(y-by);
end
end