{-# LANGUAGE LambdaCase #-}
module Stock where

import Data.Random.Normal
import System.Random
import Data.Random
import Diff
import Numeric.SpecFunctions
import Numeric.LinearAlgebra
import Util

ndigits :: Int
ndigits = 4

roundToN :: Int -> Double -> Double
roundToN n number = (fromInteger (round (number * (10^n)))) / (10.0^n)

-- constants
r :: Double
r = 0.04

variance :: Double
variance = 0.01

period :: Double
period = 360.0

yperiod :: Double
yperiod = period / 360.0

s0 :: Double
-- s0 = 50.0000
s0 = 90.000

k :: Double
-- k = 50.0000
k = 100.00

d1 :: Double
d1 = (log (s0 / k) + (r + variance / 2) * yperiod) / (sqrt (variance * yperiod))

d2 :: Double
d2 = d1 - sqrt (variance * yperiod)

-- c = s_0 * N(D_1) - exp (-r * T) * K * N(D_2)
n :: Double -> Double
n = cdf StdNormal

c :: Double
c = s0 * (n d1) - (exp (-r * yperiod)) * k * (n d2)

p :: Double
p = c + k * exp (-r * yperiod) - s0



--- basis functions
l :: Int -> String -> DiffLang
l n x
  | n >= 0    = (Exp ((Const (-0.5)) :*: (Var x))) :*: (Exp (Var x)) :*: (Const (1/(factorial n))) :*: (derivN n (laguerre n x) x)
  | otherwise = error "undefined basis function"


-- tdiff = t_{d+1} - t_d
tdiff :: Double
tdiff = 1.0/360

rzs :: Integer -> [Double]
rzs rand = normals' (0.0,1.0) (mkStdGen (fromIntegral rand))

sRand :: Integer -> Integer -> Double
sRand j n =
  let zs = rzs j
  in s zs n

s :: [Double] -> Integer -> R
s zs 0 = s0
s zs n
  | n > 0 = roundToN ndigits $ (s zs (n - 1)) * (exp $ (r - variance/2) * tdiff + (sqrt (variance * tdiff) * (zs!!(fromInteger n))))


cf :: [R] -> [R]
cf = map (\stdj -> max (k - stdj) 0)

cv :: R -> R
cv = (* (exp (-r/360.0)))

eval :: DiffLang -> R -> R
eval (Var _) x = x
eval (Const n) _ = n
eval (e1 :+: e2) x = eval e1 x + eval e2 x
eval (e1 :*: e2) x = eval e1 x * eval e2 x
eval (e1 :^: e2) x = eval e1 x ** eval e2 x
eval (Exp e) x     = exp (eval e x)

-- compute coefficients {beta_i} using (n-1)'th day's cash as predicator values,
-- n'th days cv as ys
regColN :: ([Vector R], [Vector R]) -> Int -> Matrix R
regColN (paths, cfs) n
  | length paths > n && n > 1 = regCol n laguerreBasis s t
  | otherwise = error "too big a N"
    where st :: ([R], [R], [R])
          st = unzip3 $ filter (\(a,b,c) -> a /= 0) $ zip3 (toList (cfs!!(n-1))) (toList (paths!!(n-1))) (toList (cfs!!n))
          s :: Vector R
          s = fromList $ sndOf3 st
          t :: Vector R
          t = fromList $ thdOf3 st

laguerreJ :: Int -> R -> R
laguerreJ j = eval (laguerre j "X")

laguerreBasis :: [R -> R]
laguerreBasis = map laguerreJ [0..60]

regCol :: Int -> [R -> R] -> Vector R -> Vector R -> Matrix R -- n x 1 result
regCol n basis s t = linearSolveLS a b
  where a :: Matrix R
        a = fromColumns $ map (\j -> (fromList . map (basis!!j) . toList) s) [0..n-1]
        b :: Matrix R
        b = (asColumn . fromList . map cv) (toList t)

evalPoly :: R -> [R] -> R
evalPoly x betas = foldl (\acc beta -> acc + beta * x) 0 betas

evalPolyL :: [R] -> [R] -> [R]
evalPolyL xs betas = map (\x -> evalPoly x betas) xs

cfP :: ([Vector R], [Vector R]) -> Int -> [Vector R]
cfP (paths, cfs) n
  | length paths > n && n > 1 = [newCFNP, finalCFN]
  | otherwise = error "undefined at the date!"
  where finalCFN = fromList $ map (\(s,t) -> if s == 0.0 then 0.0 else t) (zip newCFN oldCFN)
        newCFN = map (\(s,t) -> if max s t == t then s else 0.0) (zip oldCFP cfn)
        newCFNP = fromList $ map (\(s,t) -> if max s t == s then s else 0.0) (zip oldCFP cfn)
        oldCFN = toList $ cfs!!n
        oldCFP = toList $ cfs!!(n-1)
        cfn = map (\case { (0,_) -> 0 ; (_,x) -> x }) $ zip oldCFP predCF
        predCF = map (\x -> (foldl (\acc -> \case { (n,beta) -> beta * (laguerreBasis!!n) x + acc }) 0 (betas!!(n-2)))) (toList (paths!!(n-1)))
        betas :: [[(Int,R)]]
        betas = map (\n -> zip [0..n-1] ((toList . head . toColumns . regColN (paths, cfs)) n)) [2..60]

cfGenR :: (Matrix R, Matrix R) -> [Vector R]
cfGenR (paths, cfs) = foldl (\cfs' n -> take (n - 1) cfs' ++ cfP (pathsC, cfs') n ++ drop (n + 1) cfs') cfsC (reverse [2..60])
  where pathsC = toColumns paths
        cfsC = toColumns cfs

tableaux :: Integer -> Integer -> [[R]]
tableaux rand size =
  map (\r -> map (sRand r) [0..60]) ss
  where ss :: [Integer]
        ss = seeds rand size

seeds :: Integer -> Integer -> [Integer]
seeds rand length = take (fromIntegral length) $ map (\x -> (ceiling (abs x * 100000000000))) (rzs rand)

tableauxM :: Integer -> Integer -> (Matrix R, Matrix R)
tableauxM rand size = (matrix 61 paths, matrix 61 cfs)
  where paths = concat (tableaux rand size)
        cfs   = cf paths

cfGen :: (Matrix R, Matrix R) -> Matrix R
cfGen = fromColumns . tail . cfGenR

numNonZ :: Matrix R -> Int
numNonZ = sum . map (length . (filter (/= 0)) . toList) . toColumns

aOptionValue :: Matrix R -> [R]
aOptionValue m = map (\lst -> if null lst then 0.0 else (exp (-r * ((fst (head lst)) / 360.0))) * (snd (head lst)))
                 $ map (filter (\case {(n,x) -> x /= 0}) . zip [1..60] . toList) (toRows m)

eOptionValue :: Matrix R -> [R]
eOptionValue = map (* (exp (-r * yperiod))) . cf . toList . last . toColumns

cStarMin :: Integer -> (Matrix R, Matrix R) -> R
cStarMin p (cf,m) = (/dist) . sum $ map (\i -> (vs!!(i-1) - vbar) * (qs!!(i-1) - qbar))  [1..(fromIntegral p)]
  where vs = aOptionValue cf
        qs = eOptionValue m
        vbar = sum vs / (fromIntegral (length vs))
        qbar = sum qs / (fromIntegral (length qs))
        dist = sum $ map (\i -> (qs!!(i-1) - qbar)**2) [1..(fromIntegral p)]

rbar :: Integer -> (Matrix R, Matrix R) -> R -> R
rbar pass (cf,m) c = (/(fromIntegral pass)) . sum $ map (\i -> (vs!!(i-1)) - c * ((qs!!(i-1)) - p)) [1..(fromIntegral pass)]
  where vs = aOptionValue cf
        qs = eOptionValue m

output :: Integer -> Integer -> R
output rand size = rbar size (cfs, fst m) (cStarMin size (cfs, fst m))
  where cfs = cfGen m
        m = tableauxM rand size
