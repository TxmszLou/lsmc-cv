module Stock where

import Data.Random.Normal
import System.Random
import Data.Random
import Diff
import Numeric.SpecFunctions
import Data.Matrix

ndigits :: Int
ndigits = 4

roundToN :: Int -> Double -> Double
roundToN n number = (fromInteger (round (number * (10^n)))) / (10.0^n)

-- constants
r :: Double
r = 0.05

variance :: Double
variance = 0.2

period :: Double
period = 60.0

yperiod :: Double
yperiod = period / 360.0

s0 :: Double
s0 = 50.0000

k :: Double
k = 50.0000

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

s :: [Double] -> Integer -> Double
s zs 0 = 50.0000
s zs n
  | n > 0 = roundToN ndigits $ (s zs (n - 1)) * (exp $ (r - variance/2) * tdiff + (sqrt (variance * tdiff) * (zs!!(fromInteger n))))


-- average :: Int -> Int -> Double
-- average r n = ((sum . (take n)) (rzs r)) / (fromIntegral n)

seeds :: Integer -> Integer -> [Integer]
seeds rand length = take (fromIntegral length) $ map (\x -> (ceiling (abs x * 100000000000))) (rzs rand)

-- tableauxM :: Integer -> Integer -> Matrix Double
-- tableauxM rand size = matrix size 61 (\(i,j) -> )

tableaux :: Integer -> Integer -> [[Double]]
tableaux rand size =
  map (\r -> map (sRand r) [1..60]) ss
  where ss :: [Integer]
        ss = seeds rand size
