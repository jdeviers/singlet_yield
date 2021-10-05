#!/usr/bin/awk -f

NR==1{
  a=$3;
  a_init=$3;
}
{
  print "a = ",a," $3 = ",$3,"and $3-a = ",$3-a," and ($3-a)/a_init = ", ($3-a)/a_init;
  a=$3
}
END{
  print "We should sample d1*d2=a_init at every step, but the last quantity shows we sample N*(d1*d2) at every cycle (every new generation of indices), where N is the number of cycles.";
  print "This is the cause of the upward curve and possibly where the scaling error occurs. This is the case on the number of counts and also happens on the cumulative sum."
}
