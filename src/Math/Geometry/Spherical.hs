{-# LANGUAGE RankNTypes #-}

module Math.Geometry.Spherical
    ( radians
    , toRadians
    , project
    , areaTriangle
    , albers
    , perimeterPolygon
    , littow
    , craig
    , washingtonDC
    , mecca
    , winkel3
    , mercator
    , compactness1
    , relativeCompactness
    , areaConvex
    , totalPerimeter
    ) where

import           Control.Composition
import           Data.Foldable       (Foldable)

-- | Convert a `a` from degrees to radians.
radians :: Floating a => a -> a
radians = (*(pi/180))

-- | Convert both coördinates to radians.
toRadians :: Floating a => (a, a) -> (a, a)
toRadians = both radians

sinc :: Floating a => a -> a
sinc x = sin x / x

-- | For use as a reference point
washingtonDC :: (Floating a) => (a, a)
washingtonDC = toRadians (38.9072, -77.0369)

-- | For use as a reference point
mecca :: (Floating a) => (a, a)
mecca = toRadians (21.3891, 39.8579)

-- | Littow retroazimuthal + conformal projection
littow :: Floating a => (a, a) -> (a, a)
littow (long, lat) = (sin (long - referenceLong) /cos lat, cos (long - referenceLong) * tan lat)
    where referenceLong = radians 0

-- | Craig retroazimuthal projection
craig :: (Floating a, Eq a)
      => (a, a) -- ^ Reference point given in radians
      -> ((a, a) -> (a, a))
craig referencePoint (long, lat) = (long - referenceLong, y)
    where (referenceLong, referenceLat) = referencePoint
          expr = sin lat * cos (long - referenceLong) - tan referenceLat * cos lat
          y | long - referenceLong == 0 = expr
            | otherwise = (long - referenceLong) / sin (long - referenceLong) * expr

-- | Winkel Tripel projection
winkel3 :: Floating a => (a, a) -> (a, a)
winkel3 (long, lat) = ((lambda * cos phi1 + (2 * cos lat * sin (lambda/2)/sinc alpha))/2, (lat + sin lat/sinc alpha)/2)
    where lambda = long - lambda0
          phi1 = acos $ 2 / pi
          alpha = acos $ cos lat * cos (lambda/2)
          lambda0 = radians (-77.0369)

-- | Mercator projection.
mercator :: Floating a => (a, a) -> (a, a)
mercator (long, lat) = (long - meridian, asinh (tan lat))
    where meridian = radians (-98.5795)

-- | Bonne projection.
bonne :: Floating a
      => a -- ^ Standard Parallel. If you are unsure of what to put, try @radians 45@.
      -> a -- ^ Central meridian. If you are unsure of what to put, try @snd washingtonDC@.
      -> (a, a)
      -> (a, a)
bonne phi1 meridian (long, lat) = (rho * sin e, cot phi1 - rho * cos e)
    where rho = cot phi1 + phi1 - lat
          e = (long - meridian) * cos lat / rho
          cot = (1/) . tan

-- | Albers projection for a given reference point.
--
-- > ablers washingtonDC
albers :: Floating a => (a, a) -- ^ A reference point no the sphere
                     -> ((a, a) -> (a, a))
albers referencePoint (long, lat) = (rho * sin theta, rho' - rho * cos theta)
    where n = (sin phi1 + sin phi2) / 2
          theta = n * (long - referenceLong)
          c = cos phi1^(2 :: Int) + 2 * n * sin phi1
          rho = sqrt (c - 2 * n * sin lat) / n
          rho' = sqrt (c - 2 * n * sin referenceLat) / n
          -- standard parallels @ 20, 50 degrees
          phi1 = radians 20
          phi2 = radians 50
          (referenceLong, referenceLat) = referencePoint

-- | Project a given polygon
project :: (Floating a, Functor f) => ((a, a) -> (a, a)) -- ^ A projection
                                   -> f (a, a) -- ^ A polygon defined by points on the sphere
                                   -> f (a, a)
project f = fmap (f . toRadians)

-- | Compute the area of a triangle using L\'Huillier\'s formula
areaTriangle :: Floating a => a -- ^ Radius of the sphere
                           -> (a, a) -- ^ A point given in radians
                           -> (a, a)
                           -> (a, a)
                           -> a
areaTriangle r x1 x2 x3 = r^(2 :: Int) * e
    where e = 4 * atan(sqrt(tan(s/2) * tan((s - a)/2) * tan((s - b)/2) * tan((s - c)/2)))
          s = (a + b + c) / 2
          a = distanceRad x1 x2
          b = distanceRad x1 x3
          c = distanceRad x2 x3
          distanceRad = on centralAngle toRadians

-- TODO mandelbrot/fractal dimension?
-- consider "area of largest circumscribable circle" as well.

-- | Relative compactness. Dimensionless.
relativeCompactness :: (Floating a) => [(a, a)] -> a
relativeCompactness = (*scale) . compactness1
    where scale = 1/4*pi

-- | Take the area of the polygon and divide by the perimeter squared. Dimensionless.
compactness1 :: (Floating a) => [(a, a)] -> a
compactness1 p = areaPolygon p/perimeterPolygon p^(2 :: Int)

-- | Compute the area of a convex polygon on the surface of a sphere.
areaConvex :: (Floating a)
           => a -- ^ Radius of a sphere
           -> [(a, a)] -- ^ Polygon on the surface of the sphere
           -> a
areaConvex r (base1:base2:pts) = fst $ foldr stepArea (0,base2) pts
    where stepArea point (sum', base) = (sum' + areaTriangle r base1 base point, point)
areaConvex _ _ = error "attempted to take area of polygon with < 3 points"

-- | Uses areal projection; then finds area of the polygon.
-- Result is in km^2
areaPolygon :: Floating a => [(a, a)] -> a
areaPolygon = (*factor) . areaPolyRectangular . fmap (bonne (radians 25) (snd washingtonDC) . toRadians)
    where factor = 1717856/4.219690791828533e-2

-- | Given a list of polygons, return the total perimeter.
totalPerimeter :: (Functor t, Foldable t, Floating a) => t [(a, a)] -> a
totalPerimeter = sum . fmap perimeterPolygon

perimeterPolygon :: Floating a => [(a, a)] -> a
perimeterPolygon [x1, x2]       = distance x1 x2
perimeterPolygon (x1:x2:points) = perimeterPolygon (x2:points) + distance x1 x2
perimeterPolygon _              = error "Attempted to take perimeter of polygon with no points"

-- | Find the area of a polygon with rectangular coördinates given. This
-- function will throw an error when a polygon with no points is given.
areaPolyRectangular :: (Floating a) => [(a, a)] -> a
areaPolyRectangular (pt:pts) = abs . (*0.5) . fst $ foldr areaPolyCalc (0,pt) pts
    where areaPolyCalc (x2, y2) (sum',(x1,y1)) = (sum' + (x1 * y2 - x2 * y1),(x2,y2))
areaPolyRectangular _ = error "Attempted to take area of polygon with no points"

-- | Distance in kilometers between two points given in degrees.
distance :: Floating a => (a, a) -> (a, a) -> a
distance = (*6371) .* on centralAngle toRadians

-- | Compute central angle from points given in radians
centralAngle :: Floating a => (a, a) -> (a, a) -> a
centralAngle (long1, lat1) (long2, lat2) =
    acos $ sin lat1 * sin lat2 + cos lat1 * cos lat2 * cos (long1 - long2)
