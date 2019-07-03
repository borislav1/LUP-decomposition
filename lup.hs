import Data.List 

outer_product v wT = map (\x -> map (x*) wT) v
subtract_vectors v w = zipWith (-) v w
subtract_matrices m1 m2 = zipWith subtract_vectors m1 m2

schur_complement mat = (schur_mat, v)
                where x = head (head mat)
                      wT = tail (head mat)
                      v = map (\ y -> y / x) $ reverse $ foldl (\ acc xs -> (head xs):acc) [] (tail mat)
                      mat' = reverse $ foldl (\ acc xs -> (tail xs):acc) [] (tail mat)
                      schur_mat = subtract_matrices mat' (outer_product v wT)

                    
max_pivot mat = max_pivot' mat (head $ head mat) 0 0 
max_pivot' [] _ _ index = index
max_pivot' (v:mat) max counter index | (abs $ head v) > max = max_pivot' mat (abs $ head v) (counter+1) counter 
                                     | otherwise = max_pivot' mat max (counter+1) index


identity_matrix n = [ xs ++ [1] ++ ys | x <- ns, y <- (reverse ns), let xs = f x, let ys = f y, x + y == (n-1)]
        where f m = replicate m 0
              ns = [0..(n-1)] 

make_swap_perm i j n = swapElementsAt i j (identity_matrix n)

mm_mul a b = [ [ sum $ zipWith (*) ar bc | bc <- (transpose b) ] | ar <- a ]
mv_mul a b = concat $ mm_mul a (map (\x -> [x]) b)


lup_decompose mat = lup_decompose' mat

lup_decompose' [x] = ([[1]], [x], [[1]])
lup_decompose' a = ((zs:(zipWith (:) (mv_mul p' v) l')), ((head a'):(map (\xs -> 0:xs) u')), ((embed p') `mm_mul` q) )
             where (l', u', p') = lup_decompose' mat'                   
                   pivot = max_pivot a
                   q = make_swap_perm 0 pivot (length a)
                   a' = q `mm_mul` a
                   (mat', v) = schur_complement a'
                   zs = 1 : (replicate (length mat') 0)
                   embed per = zs:(map (\xs -> 0:xs) per) 


forward_substitution (xs:lower) (b:bs) = forward_substitution' lower bs 2 [b / head xs]

forward_substitution' [] _ _ ys = ys  
forward_substitution' (xs:lower) (b:bs) n ys = forward_substitution' lower bs (n+1) (ys++[x]) 
                                     where v = take (length v'-1) v'
                                           v' = take n xs  
                                           y = last v'
                                           x = (-(sum $ zipWith (*) v ys) + b) / y

back_substitution upper ys = reverse $ forward_substitution (reverse $ map reverse upper) (reverse ys)  

solve_equation mat b = back_substitution upper y
                    where (lower, upper, p) = lup_decompose mat 
                          y = forward_substitution lower (p `mv_mul` b)


swapElementsAt i j xs | i == j = xs
                      | otherwise  = left ++ [elemJ] ++ middle ++ [elemI] ++ right
                        where elemI = xs !! i
                              elemJ = xs !! j
                              left = take i xs
                              middle = take (j - i - 1) (drop (i + 1) xs)
                              right = drop (j + 1) xs

          
