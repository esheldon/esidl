FUNCTION pq_limits, x

  return, [-sqrt(!conv_rmax^2-x^2),sqrt(!conv_rmax^2-x^2)]
END 
