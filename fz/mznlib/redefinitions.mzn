% Mistral redefinitions 

predicate int_eq_imp(var int: a, var int: b, var bool: r);
predicate int_ne_imp(var int: a, var int: b, var bool: r);
predicate int_ge_imp(var int: a, var int: b, var bool: r);
predicate int_gt_imp(var int: a, var int: b, var bool: r);
predicate int_le_imp(var int: a, var int: b, var bool: r);
predicate int_lt_imp(var int: a, var int: b, var bool: r);

predicate count_eq_reif(array[int] of var int: x, var int: y, var int: c, var bool: b);
predicate count_leq_reif(array[int] of var int: x, var int: y, var int: c, var bool: b);

predicate array_int_maximum(var int: m, array[int] of var int: x);
predicate array_int_minimum(var int: m, array[int] of var int: x);

