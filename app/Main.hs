module Main where

import Stock

main :: IO ()
main =
  (putStrLn . show) $ tableaux 24725231 50000
