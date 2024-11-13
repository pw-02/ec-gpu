// Elliptic curve operations (Short Weierstrass Jacobian form)

#define POINT_ZERO ((POINT_projective){FIELD_ZERO, FIELD_ONE, FIELD_ZERO})

typedef struct {
  FIELD x;
  FIELD y;
} POINT_affine;

typedef struct {
  FIELD x;
  FIELD y;
  FIELD z;
} POINT_projective;


// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
DEVICE POINT_projective POINT_double(POINT_projective inp) {
  //const FIELD local_zero = FIELD_ZERO;
  //if(FIELD_eq(inp.z, local_zero)) {
  //   return inp;
  //}

  const FIELD local_zero = FIELD_ZERO;
  const FIELD local_one = FIELD_ONE;
  const FIELD local_two = FIELD_add(local_one,local_one);
  const FIELD local_three = FIELD_add(local_two,local_one);
  const FIELD local_six = FIELD_add(local_three,local_three);
  const FIELD local_nine = FIELD_add(local_six,local_three);


  FIELD t0 = FIELD_sqr(inp.y); // let t0 = self.y.square();
  FIELD z3 =FIELD_add(t0, t0); // let z3 = t0 + t0;
   z3 = FIELD_add(z3, z3); //let z3 = z3 + z3;
   z3 = FIELD_add(z3, z3); //let z3 = z3 + z3;
  FIELD t1 = FIELD_mul(inp.y,inp.z);  //let t1 = self.y * self.z;
  FIELD t2 = FIELD_sqr(inp.z); //let t2 = self.z.square();
  t2 = FIELD_mul(t2, local_nine);  //let t2 = $name::mul_by_3b(&t2);

  //t2 = FIELD_mul(t2, t2);  //let t2 = $name::mul_by_3b(&t2);

  FIELD x3 = FIELD_mul(t2,z3); //let x3 = t2 * z3;
  FIELD y3 = FIELD_add(t0,t2); // let y3 = t0 + t2;
   z3 = FIELD_mul(t1,z3); // let z3 = t1 * z3;
   t1 = FIELD_add(t2,t2); // let t1 = t2 + t2;
   t2 = FIELD_add(t1,t2); // let t2 = t1 + t2;
   t0 = FIELD_sub(t0,t2); //  let t0 = t0 - t2;
   y3 = FIELD_mul(t0,y3); // let y3 = t0 * y3;
   y3  = FIELD_add(x3,y3); //  let y3 = x3 + y3;
   t1 =  FIELD_mul(inp.x,inp.y); //let t2 = self.z.square();
   x3 = FIELD_mul(t0,t1); // let x3 = t0 * t1;
   x3 = FIELD_add(x3,x3);

  POINT_projective tmp;
  tmp.x = x3;
  tmp.y = y3;
  tmp.z = z3;
  return tmp;

  if(FIELD_eq(inp.z, local_zero)) {
      return inp;
  }else
  {
  return tmp;
  }
 
}


/*
// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
DEVICE POINT_projective POINT_double(POINT_projective inp) {
  const FIELD local_zero = FIELD_ZERO;
  if(FIELD_eq(inp.z, local_zero)) {
      return inp;
  }

  const FIELD a = FIELD_sqr(inp.x); // A = X1^2
  const FIELD b = FIELD_sqr(inp.y); // B = Y1^2
  FIELD c = FIELD_sqr(b); // C = B^2

  // D = 2*((X1+B)2-A-C)
  FIELD d = FIELD_add(inp.x, b);
  d = FIELD_sqr(d); d = FIELD_sub(FIELD_sub(d, a), c); d = FIELD_double(d);

  const FIELD e = FIELD_add(FIELD_double(a), a); // E = 3*A
  const FIELD f = FIELD_sqr(e);

  inp.z = FIELD_mul(inp.y, inp.z); inp.z = FIELD_double(inp.z); // Z3 = 2*Y1*Z1
  inp.x = FIELD_sub(FIELD_sub(f, d), d); // X3 = F-2*D

  // Y3 = E*(D-X3)-8*C
  c = FIELD_double(c); c = FIELD_double(c); c = FIELD_double(c);
  inp.y = FIELD_sub(FIELD_mul(FIELD_sub(d, inp.x), e), c);

  return inp;
}
*/


DEVICE POINT_projective POINT_add_mixed(POINT_projective a, POINT_affine b) {
  
  const FIELD local_zero = FIELD_ZERO;
  const FIELD local_one = FIELD_ONE;
  const FIELD local_two = FIELD_add(local_one,local_one);
  const FIELD local_three = FIELD_add(local_two,local_one);
  const FIELD local_six = FIELD_add(local_three,local_three);
  const FIELD local_nine = FIELD_add(local_six,local_three);

  FIELD t0 = FIELD_mul(a.x,b.x);  //let t0 = self.x * rhs.x;
  FIELD t1 = FIELD_mul(a.y,b.y); // let t1 = self.y * rhs.y;
  FIELD t3 = FIELD_add(b.x,b.y); // let t3 = rhs.x + rhs.y;
  FIELD t4 = FIELD_add(a.x,a.y); // let t4 = self.x + self.y;
  t3 = FIELD_mul(t3,t4); // let t3 = t3 * t4;
  t4 = FIELD_add(t0,t1); // let t4 = t0 + t1;
  t3 = FIELD_sub(t3,t4); // let t3 = t3 - t4;
  t4 = FIELD_mul(b.y,a.z); // let t4 = rhs.y * self.z;
  t4 = FIELD_add(t4,a.y); // let t4 = t4 + self.y;
  FIELD y3 = FIELD_mul(b.x,a.z); // let y3 = rhs.x * self.z;
  y3 = FIELD_add(y3,a.x); // let y3 = y3 + self.x;
  FIELD x3 = FIELD_add(t0,t0); // let x3 = t0 + t0;
  t0 = FIELD_add(x3,t0); // let t0 = x3 + t0;
  FIELD t2 = FIELD_mul(a.z, local_nine); // let t2 = $name::mul_by_3b(&self.z);
  FIELD z3 = FIELD_add(t1,t2); // let z3 = t1 + t2;
  t1 = FIELD_sub(t1,t2); // let t1 = t1 - t2;
  y3 = FIELD_mul(y3,local_nine); // let y3 = $name::mul_by_3b(&y3);
  x3 = FIELD_mul(t4,y3); // let x3 = t4 * y3;
  t2 = FIELD_mul(t3,t1); // let t2 = t3 * t1;
  x3 = FIELD_sub(t2,x3); // let x3 = t2 - x3;
  y3 = FIELD_mul(y3,t0); // let y3 = y3 * t0;
  t1 = FIELD_mul(t1,z3); // let t1 = t1 * z3;
  y3 = FIELD_add(t1,y3); // let y3 = t1 + y3;
  t0 = FIELD_mul(t0,t3); // let t0 = t0 * t3;
  z3 = FIELD_mul(z3,t4); // let z3 = z3 * t4;
  z3 = FIELD_add(z3,t0); // let z3 = z3 + t0;

   POINT_projective tmp;
   tmp.x = x3;
   tmp.y = y3;
   tmp.z = z3;
  
  //$name::conditional_select(&tmp, self, rhs.is_identity())

    if(FIELD_eq(b.x, local_zero)) 
    {
      if(FIELD_eq(b.y, local_zero))
      {
        return a;
      }
      else
      {
        return tmp;
      }
    } else
    {
       return tmp;
    }
}

/*
// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
DEVICE POINT_projective POINT_add_mixed(POINT_projective a, POINT_affine b) {
  const FIELD local_zero = FIELD_ZERO;
  if(FIELD_eq(a.z, local_zero)) {
    const FIELD local_one = FIELD_ONE;
    a.x = b.x;
    a.y = b.y;
    a.z = local_one;
    return a;
  }

  const FIELD z1z1 = FIELD_sqr(a.z);
  const FIELD u2 = FIELD_mul(b.x, z1z1);
  const FIELD s2 = FIELD_mul(FIELD_mul(b.y, a.z), z1z1);

  if(FIELD_eq(a.x, u2) && FIELD_eq(a.y, s2)) {
      return POINT_double(a);
  }

  const FIELD h = FIELD_sub(u2, a.x); // H = U2-X1
  const FIELD hh = FIELD_sqr(h); // HH = H^2
  FIELD i = FIELD_double(hh); i = FIELD_double(i); // I = 4*HH
  FIELD j = FIELD_mul(h, i); // J = H*I
  FIELD r = FIELD_sub(s2, a.y); r = FIELD_double(r); // r = 2*(S2-Y1)
  const FIELD v = FIELD_mul(a.x, i);

  POINT_projective ret;

  // X3 = r^2 - J - 2*V
  ret.x = FIELD_sub(FIELD_sub(FIELD_sqr(r), j), FIELD_double(v));

  // Y3 = r*(V-X3)-2*Y1*J
  j = FIELD_mul(a.y, j); j = FIELD_double(j);
  ret.y = FIELD_sub(FIELD_mul(FIELD_sub(v, ret.x), r), j);

  // Z3 = (Z1+H)^2-Z1Z1-HH
  ret.z = FIELD_add(a.z, h); ret.z = FIELD_sub(FIELD_sub(FIELD_sqr(ret.z), z1z1), hh);
  return ret;
}
*/

// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
DEVICE POINT_projective POINT_add(POINT_projective a, POINT_projective b) {

  const FIELD local_zero = FIELD_ZERO;
  const FIELD local_one = FIELD_ONE;
  const FIELD local_two = FIELD_add(local_one,local_one);
  const FIELD local_three = FIELD_add(local_two,local_one);
  const FIELD local_six = FIELD_add(local_three,local_three);
  const FIELD local_nine = FIELD_add(local_six,local_three);

  FIELD t0 = FIELD_mul(a.x, b.x);//let t0 = self.x * rhs.x;
  FIELD t1 = FIELD_mul(a.y, b.y);	//let t1 = self.y * rhs.y;
  FIELD t2 = FIELD_mul(a.z, b.z);	//let t2 = self.z * rhs.z;
  FIELD t3 = FIELD_add(a.x, a.y); //let t3 = self.x + self.y;
  FIELD t4 = FIELD_add(b.x,b.y); //let t4 = rhs.x + rhs.y;
   t3 = FIELD_mul(t3,t4);//let t3 = t3 * t4;
   t4 = FIELD_add(t0,t1);//let t4 = t0 + t1;
   t3 = FIELD_sub(t3,t4);//let t3 = t3 - t4;
   t4 = FIELD_add(a.y,a.z);//let t4 = self.y + self.z;
  FIELD x3 = FIELD_add(b.y,b.z);//let x3 = rhs.y + rhs.z;
   t4 = FIELD_mul(t4,x3);//let t4 = t4 * x3;
   x3 = FIELD_add(t1,t2);//let x3 = t1 + t2;
   t4 = FIELD_sub(t4,x3);//let t4 = t4 - x3;
   x3 = FIELD_add(a.x,a.z);//let x3 = self.x + self.z;
  FIELD y3 = FIELD_add(b.x,b.z);//let y3 = rhs.x + rhs.z;
   x3 = FIELD_mul(x3,y3);//let x3 = x3 * y3;
   y3 = FIELD_add(t0,t2);//let y3 = t0 + t2;
   y3 = FIELD_sub(x3,y3);//let y3 = x3 - y3;
   x3 = FIELD_add(t0,t0);//let x3 = t0 + t0;
   t0 = FIELD_add(x3,t0);//let t0 = x3 + t0;
   t2 = FIELD_mul(t2, local_nine); // let t2 = $name::mul_by_3b(&t2);
   FIELD z3 = FIELD_add(t1,t2);//let z3 = t1 + t2;
   t1 = FIELD_sub(t1,t2);//let t1 = t1 - t2;
   y3 = FIELD_mul(y3, local_nine); // let y3 = $name::mul_by_3b(&y3);
   x3 = FIELD_mul(t4,y3);//let x3 = t4 * y3;
   t2 = FIELD_mul(t3,t1);//let t2 = t3 * t1;
   x3  = FIELD_sub(t2,x3);//let x3 = t2 - x3;
   y3 = FIELD_mul(y3,t0);//let y3 = y3 * t0;
   t1 = FIELD_mul(t1,z3);//let t1 = t1 * z3;
   y3  = FIELD_add(t1,y3);//let y3 = t1 + y3;
   t0 = FIELD_mul(t0,t3);//let t0 = t0 * t3;
   z3 = FIELD_mul(z3,t4);//let z3 = z3 * t4;
   z3= FIELD_add(z3,t0);//let z3 = z3 + t0;

   POINT_projective tmp;
   tmp.x = x3;
   tmp.y = y3;
   tmp.z = z3;
   return tmp;
  
}

/*
// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
DEVICE POINT_projective POINT_add(POINT_projective a, POINT_projective b) {

  const FIELD local_zero = FIELD_ZERO;
  if(FIELD_eq(a.z, local_zero)) return b;
  if(FIELD_eq(b.z, local_zero)) return a;

  const FIELD z1z1 = FIELD_sqr(a.z); // Z1Z1 = Z1^2
  const FIELD z2z2 = FIELD_sqr(b.z); // Z2Z2 = Z2^2
  const FIELD u1 = FIELD_mul(a.x, z2z2); // U1 = X1*Z2Z2
  const FIELD u2 = FIELD_mul(b.x, z1z1); // U2 = X2*Z1Z1
  FIELD s1 = FIELD_mul(FIELD_mul(a.y, b.z), z2z2); // S1 = Y1*Z2*Z2Z2
  const FIELD s2 = FIELD_mul(FIELD_mul(b.y, a.z), z1z1); // S2 = Y2*Z1*Z1Z1

  if(FIELD_eq(u1, u2) && FIELD_eq(s1, s2))
    return POINT_double(a);
  else {
    const FIELD h = FIELD_sub(u2, u1); // H = U2-U1
    FIELD i = FIELD_double(h); i = FIELD_sqr(i); // I = (2*H)^2
    const FIELD j = FIELD_mul(h, i); // J = H*I
    FIELD r = FIELD_sub(s2, s1); r = FIELD_double(r); // r = 2*(S2-S1)
    const FIELD v = FIELD_mul(u1, i); // V = U1*I
    a.x = FIELD_sub(FIELD_sub(FIELD_sub(FIELD_sqr(r), j), v), v); // X3 = r^2 - J - 2*V

    // Y3 = r*(V - X3) - 2*S1*J
    a.y = FIELD_mul(FIELD_sub(v, a.x), r);
    s1 = FIELD_mul(s1, j); s1 = FIELD_double(s1); // S1 = S1 * J * 2
    a.y = FIELD_sub(a.y, s1);

    // Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
    a.z = FIELD_add(a.z, b.z); a.z = FIELD_sqr(a.z);
    a.z = FIELD_sub(FIELD_sub(a.z, z1z1), z2z2);
    a.z = FIELD_mul(a.z, h);

    return a;
  }
  
} */

