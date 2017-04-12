module Game where

import Matrix

-- nxn active strategy matrices only
game_SLE :: [[Double]] -> (Double, [Double], [Double])
game_SLE a = (nu, map (nu*) x, map (nu*) y) where  v = [1|_ <- [1..(length a)]]
                                                   y = cramer a v
                                                   x = cramer (transpose a) v
                                                   nu = (sum x) ** (-1)
-- nxn active strategy matrices only
game_INV :: [[Double]] -> (Double, [Double], [Double])
game_INV a = (nu, map (nu*) $ head x', map (nu*) $ head y') where a' = inv a
                                                                  x' = mtxMult [[1|_ <- [1..(length a)]]] a'
                                                                  y' = transpose $ mtxMult a' [[1]|_ <- [1..(length a)]]
                                                                  nu = (/) 1 $ sum $ head y'

activeStrat :: Ord a => [[a]] -> [[a]]
activeStrat a =  if length a /= length filtered then transpose $ activeStrat $ transpose filtered else a where leq x y = forAll id $ zipWith (<=) x y
                                                                                                               filtered = [head $ drop n a| n <- [0..(length a - 1)], not $ exists (leq $ head $ drop n a) $ take n a ++ drop (n + 1) a]
