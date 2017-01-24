module Main where

import Stock

main :: IO ()
main = do
  rand <- getLine
  size <- getLine
  putStrLn $ show $ output (read rand) (read size)
