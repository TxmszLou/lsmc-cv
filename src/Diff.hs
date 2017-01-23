module Diff where

data DiffLang
  = Var String
  | Const Double
  | DiffLang :+: DiffLang
  | DiffLang :*: DiffLang
  | DiffLang :^: DiffLang
  | Exp DiffLang
  deriving Eq

instance Show DiffLang where
  show (Var x) = x
  show (Const n) = show n
  show (e1 :+: e2) = "(" ++ show e1 ++ "*" ++ show e2 ++ ")"
  show (e1 :*: e2) = "(" ++ show e1 ++ "*" ++ show e2 ++ ")"
  show (e1 :^: e2) = "(" ++ show e1 ++ "^" ++ show e2 ++ ")"
  show (Exp e)     = "(e^" ++ show e ++ ")"

laguerre :: Int -> String -> DiffLang
laguerre n x = ((Var x) :^: (Const (fromIntegral n))) :*: (Exp ((Const (-1)) :*: (Var x)))

derivN :: Int -> DiffLang -> String -> DiffLang
derivN 0 e x = e
derivN n e x
  | n > 0     = derivN (n - 1) (deriv e x) x
  | otherwise = error "differentiate negative times"

deriv :: DiffLang -> String -> DiffLang
deriv e x = simplify (deriv_base e x)

deriv_base :: DiffLang -> String -> DiffLang
deriv_base (Var x) x'
  | x == x'   = Const 1
  | otherwise = Const 1
deriv_base (Const n) _ = Const 0
deriv_base (e1 :+: e2) x = (deriv_base e1 x) :+: (deriv_base e2 x)
deriv_base (e1 :*: e2) x = ((deriv_base e1 x) :*: e2) :+: (e1 :*: (deriv_base e2 x))
deriv_base (e1 :^: e2) x = ((Const n) :*: (e1 :^: (Const (n - 1)))) :*: (deriv_base e1 x)
  where n =
          case simplify e2 of
            Const m -> m
            _       -> error "differentiating wrt non constant!"
deriv_base (Exp e) x = (Exp e) :*: (deriv_base e x)

simplify :: DiffLang -> DiffLang
simplify (Var x) = Var x
simplify (Const n) = Const n
simplify (e1 :+: e2) = simplify_sum (simplify e1 :+: simplify e2)
simplify (e1 :*: e2) = simplify_prod (simplify e1 :*: simplify e2)
simplify (e1 :^: e2) = simplify_pow (simplify e1 :^: simplify e2)
simplify (Exp e)     = simplify_exp (Exp (simplify e))

simplify_sum :: DiffLang -> DiffLang
simplify_sum ((Const 0.0) :+: e) = e
simplify_sum (e :+: (Const 0.0)) = e
simplify_sum ((Const n) :+: (Const m)) = Const (n + m)
simplify_sum (e1 :+: e2) = (e1 :+: e2)

simplify_prod :: DiffLang -> DiffLang
simplify_prod ((Const 0.0) :*: e) = Const 0.0
simplify_prod (e :*: (Const 0.0)) = Const 0.0
simplify_prod ((Const 1.0) :*: e) = e
simplify_prod (e :*: (Const 1.0)) = e
simplify_prod ((Const n) :*: (Const m)) = Const (n * m)
simplify_prod (e1 :*: e2) = (e1 :*: e2)

simplify_pow :: DiffLang -> DiffLang
simplify_pow (_ :^: (Const 0.0)) = Const 1.0
simplify_pow (e :^: (Const 1.0)) = e
simplify_pow ((Const n) :^: (Const m))
  | ceiling m == floor m = Const (n ^^ (round m))
  | otherwise = error "taking power of non-integer"
simplify_pow (e1 :^: e2) = e1 :^: e2

simplify_exp :: DiffLang -> DiffLang
simplify_exp (Exp (Const n))   = Const (exp n)
simplify_exp e                 = e
