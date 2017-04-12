module Matrix where

detLaplace :: [[Integer]] -> Integer
detLaplace [[a]] = a
detLaplace a     = sum $ map (\x -> (head $ drop (x - 1) $ last a) * (detLaplace [take (x - 1) y ++ drop x y| y <- init a]) * ((-1)^(x + n))) [1..n] where n = length a


permute    :: [a] -> [[a]]
permute [a]      = [[a]]
permute (x:xs)   = concat $ map (\p -> [take p y ++ [x] ++ drop p y|y <- permute xs]) [0..n] where n = length xs

parity     :: [Int] -> Int
parity [a]       = 1
parity (x:xs)    = parity xs * product (map (\p -> signum $ p - x) xs)

detPerm    :: [[Int]] -> Int
detPerm   a      = sum $ map (\p -> parity p * (product $ zipWith (\c d -> head $ drop (c - 1) d) p a)) $ permute [1..n] where n = length a


transpose  :: [[a]] -> [[a]]
transpose ([]:_) = []
transpose a      = (map head a):(transpose $ map tail a)

forAll     :: (a -> Bool) -> [a] -> Bool
forAll _ []      = True
forAll f (x:xs)  = f x && forAll f xs

exists     :: (a -> Bool) -> [a] -> Bool 
exists _ []      = False
exists f (x:xs)  = f x || exists f xs

detGauss     :: (RealFloat a) => [[a]] -> a
detGauss [[a]]   = a
detGauss a@((x:_):_) | x == 0    = if exists (0/=) $ head $ transpose a then (*) (-1) $ detGauss $ head [(head $ drop y a):(take y a ++ drop (y + 1) a)|y <- [1..], head (head $ drop y a) /= 0] else 0 
                     | otherwise = (*) x $ detGauss $ map (\k -> tail $ zipWith (\q p -> p - ((q * head k)/x)) (head a) k) $ tail a

mtxMult      :: (RealFloat a) => [[a]] -> [[a]] -> [[a]]
mtxMult a b      = map (\x -> map (\y -> sum $ zipWith (*) x y) $ transpose b) a

mtxMultC     :: (RealFloat a) => a -> [[a]] -> [[a]]
mtxMultC a b     = map (map (a*)) b 

mtxAdd       :: (RealFloat a) => [[a]] -> [[a]] -> [[a]]
mtxAdd a b       = zipWith (zipWith (+)) a b

mtxSub       :: (RealFloat a) => [[a]] -> [[a]] -> [[a]]
mtxSub a b       = zipWith (zipWith (-)) a b

mtxPow       :: (RealFloat a) => [[a]] -> Int -> [[a]]
mtxPow a 0       = mtxOne $ length a
mtxPow a n       = mtxMult a $ mtxPow a $ n - 1

vcrNorm2     :: (RealFloat a) => [a] -> a
vcrNorm2 a       = sqrt $ sum $ map (**2) a

mtxNorm2     :: (RealFloat a) => [[a]] -> a
mtxNorm2 a       = vcrNorm2 $ concat a

mtxZero      :: (RealFloat a) => Int -> [[a]]
mtxZero n        = [[0|_<-k]|_<-k] where k = [1..n]

mtxOne       :: (RealFloat a) => Int -> [[a]]
mtxOne 1         = [[1]]
mtxOne n         = ([1] ++ [0|_<-[1..(n - 1)]]):(map ([0] ++) $ mtxOne $ n - 1)

-- Default way to calculate determinants. Change to whatever determinant computing function you want.
det              = detGauss 

-- squeezes a vector into an nxn matrix
mtx          :: [a] -> [[a]]
mtx v            = [take n $ drop (n * x) v|x <- [0..(n - 1)]] where n = floor $ sqrt $ fromIntegral $ length v

mtxMap       :: (a -> b) -> [[a]] -> [[b]]
mtxMap f a       = map (map f) a

vcrMax       :: Ord a => [a] -> a
vcrMax a         = head [x| x <- a, forAll (x>=) a]

vcrMin       :: Ord a => [a] -> a
vcrMin a         = head [x| x <- a, forAll (x<=) a]

mtxMax       :: Ord a => [[a]] -> a
mtxMax a         = vcrMax $ concat a

mtxMin       :: Ord a => [[a]] -> a
mtxMin a         = vcrMin $ concat a

mtxSum       :: (RealFloat a) =>  [[[a]]] -> [[a]]
mtxSum [a]       = a
mtxSum (x:xs)    = mtxAdd x $ mtxSum xs

mtxProd      :: (RealFloat a) =>  [[[a]]] -> [[a]]  
mtxProd [a]      = a
mtxProd (x:xs)   = mtxMult x $ mtxProd xs

mtxPoly      :: (RealFloat a) =>  [a] -> [[a]] -> [[a]]
mtxPoly [c] a    = mtxMultC c $ mtxOne $ length a
mtxPoly c a      = mtxAdd (mtxPoly [head c] a) $ mtxMult a $ mtxPoly (tail c) a

adj          :: [[Double]] -> [[Double]]
adj a            = [[let t q p = take (q - 1) p ++ drop q p in (*) ((^) (-1) $ x + y) $ det $ t y $ map (t x) a| y <- [1..n]] | x <- [1..n]] where n = length a


invAdj       :: [[Double]] -> [[Double]]
invAdj a           = mtxMap (/ det a) $ adj a

invSchultz   :: (RealFloat a) => Int -> a -> [[a]] -> [[a]]
invSchultz n err a = schultz n err a (mtxMultC ((/) 1 $ mtxNorm2 $ mtxMult a t) t) where t = transpose a
                                                                                         schultz m err a u = let psi = mtxSub (mtxOne (length a)) $ mtxMult a u in if mtxNorm2 psi < err then u else schultz m err a (mtxMult u $ mtxPoly [1|_<-[0..m]] psi) 

invGaussJordan :: (RealFloat a) => [[a]] -> [[a]]
invGaussJordan a0  = gauss $ jordan (a0, mtxOne $ length a0) where takeLast n xs = foldl (const . tail) xs (drop n xs)
{-                   ^ lol why not^       -}                       dropLast n xs = foldl (const . init) xs [1..n]
                                                                   gauss ([[1]], b) = b
                                                                   gauss (a, b)     = gauss (map init $ init a, (zipWith (\x -> zipWith (\y z -> z - y*x) current) (init $ map last a) $ take (length a - 1) b) ++ (drop (length a - 1) b)) where current = last $ take (length a) b
                                                                   jordan  ([[a]], b) = ([[1]], init b ++ [map (/a) $ last b])
                                                                   jordan  (a@((x:_):_), b) | x == 0    = jordan $ head [let swap m = (take (length m - length a) m) ++ [last $ dropLast y m] ++ (drop (length m - length a) $ dropLast (y + 1) m ++ takeLast y m) in (swap a, swap b)|y <- [0..], head (last $ dropLast y a) /= 0]
                                                                                            | otherwise = ((head a'):(map (0:) $ fst next), snd next) where (a', b')   = ((map (/x) $ head a):(tail a), dropLast n b ++ [map (/x) $ head $ takeLast n b] ++ (takeLast (n - 1) b))
                                                                                                                                                            (a'', b'') = (map (\k -> tail $ zipWith (\q p -> p - q * head k) (head a') k) $ tail a', dropLast (n - 1) b' ++ (zipWith (\k l -> zipWith (\q p -> p - q * head l) (head $ takeLast n b') k) (takeLast (n - 1) b') $ tail a'))
                                                                                                                                                            n          = length a
                                                                                                                                                            next       = jordan (a'', b'')

-- Default way to find inverse matrices. Change to whatever inversion you want.
inv                = invGaussJordan


cramer       :: [[Double]] -> [Double] -> [Double]
cramer a0 v      = [(/) (det (take (y - 1) a ++ [v] ++ drop y a)) $ det a|y <- [1..(length a)]] where a = transpose a0 


opApply      :: (RealFloat a) => [[a]] -> [a] -> [a]
opApply o v      = head $ transpose $ mtxMult o $ transpose [v]


