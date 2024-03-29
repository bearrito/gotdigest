package gotdigest

import (
	"math"
	"slices"
)

type Centroid struct {
	size int
	mean float64
}

type TDigest struct {
	centroids     []Centroid
	delta         float64
	shouldBuffer  bool
	bufferMaxSize int
	buffer        []float64
	useNormalized bool
}

func last(s *[]Centroid) *Centroid {

	return &((*s)[len(*s)-1])
}

type unNormalizedPotential func(float64, float64) float64
type normalizedPotential func(float64, float64, int) float64

func k0Potential(quantile float64, delta float64) float64 {

	scalar := delta / (2.0)
	inverse := quantile
	return scalar * inverse
}

func k1Potential(quantile float64, delta float64) float64 {

	scalar := delta / (2.0 * math.Pi)
	inverse := math.Asin(2.0*quantile - 1.0)
	return scalar * inverse
}

func k2Potential(quantile float64, delta float64, n int) float64 {

	if quantile == 0 {
		quantile = 0.000000001

	}

	if quantile == 1 {

		quantile = .999999999
	}

	denom := (4.0*math.Log(float64(n)/delta) + 24.0)
	lhs := delta / denom
	rhs := math.Log(quantile / (1.0 - quantile))

	return lhs * rhs
}

func NewCentroidWithValue(val float64) Centroid {

	return Centroid{1, val}
}

func NewDigestFromValues(vals ...float64) TDigest {
	centroids := make([]Centroid, len(vals))
	digest := TDigest{centroids: centroids}
	for val := range vals {
		digest.centroids = append(digest.centroids, NewCentroidWithValue(float64(val)))

	}

	return digest

}

func (c *Centroid) weight() float64 {

	return c.mean * float64(c.size)
}

func (c *Centroid) update(other Centroid) {

	weight := c.weight() + other.weight()
	c.size += other.size
	c.mean = weight / float64(c.size)
}

func mergeCentroids(first Centroid, second Centroid) Centroid {

	size := (first.size + second.size)
	weight := first.weight() + second.weight()

	return Centroid{size, weight / float64(size)}

}

func NewDigest(delta float64, shouldBuffer bool, bufferMax int, useNormalized bool) TDigest {

	bins := make([]Centroid, 0)
	buffer := make([]float64, 0, bufferMax)
	return TDigest{bins, delta, shouldBuffer, bufferMax, buffer, useNormalized}

}

func NewDigestFromBin(bin Centroid) TDigest {

	bins := make([]Centroid, 1)
	bins[0] = bin

	return TDigest{centroids: bins}

}

func (d *TDigest) append(value float64) {

	if d.shouldBuffer {
		d.buffer = append(d.buffer, value)
		if len(d.buffer) == d.bufferMaxSize {
			slices.Sort(d.buffer)

			centroids := make([]Centroid, 0, len(d.buffer))

			for _, data := range d.buffer {
				centroids = append(centroids, Centroid{1, data})

			}
			d.buffer = make([]float64, 0, d.bufferMaxSize)

			digest := TDigest{centroids: centroids}
			d.merge(&digest)
			d.centroids = compressBins(d.centroids, d.delta)
		}

	} else {
		b := NewCentroidWithValue(value)
		digest := NewDigestFromBin(b)
		d.merge(&digest)
		d.centroids = compressBins(d.centroids, d.delta)
	}

}

func (d *TDigest) count() int {

	count := 0
	for _, bin := range d.centroids {
		count += bin.size

	}
	return count

}

func (d *TDigest) merge(other *TDigest) {

	size := len(d.centroids) + len(other.centroids)
	merged := make([]Centroid, 0, size)

	dCount := 0
	otherCount := 0
	for dCount < len(d.centroids) && otherCount < len(other.centroids) {

		if d.centroids[dCount].mean < other.centroids[otherCount].mean {
			merged = append(merged, d.centroids[dCount])
			dCount++
		} else {

			merged = append(merged, other.centroids[otherCount])
			otherCount++

		}

	}

	for dCount < len(d.centroids) {
		merged = append(merged, d.centroids[dCount])
		dCount++
	}

	for otherCount < len(other.centroids) {
		merged = append(merged, other.centroids[otherCount])
		otherCount++
	}

	d.centroids = merged
}

func compressBins(centroids []Centroid, delta float64) []Centroid {

	if len(centroids) == 0 {
		return centroids
	}

	totalSize := 0
	for _, centroid := range centroids {
		totalSize += centroid.size

	}

	compressedCentroids := make([]Centroid, 0)
	compressedCentroids = append(compressedCentroids, centroids[0])

	accumulatedSize := compressedCentroids[0].size

	minPotential := k1Potential(0, delta)
	i := 1
	for i < len(centroids) {
		nextCentroid := centroids[i]
		i++

		quotientIndex := float64((accumulatedSize + nextCentroid.size)) / float64(totalSize)
		if quotientIndex > 1.0 {
			panic("Cannot have quantiles greater than 1.0")
		}
		if k1Potential(quotientIndex, delta)-minPotential <= 1 {
			last(&compressedCentroids).update(nextCentroid)
		} else {
			compressedCentroids = append(compressedCentroids, nextCentroid)
			quantile := float64(accumulatedSize) / float64(totalSize)
			minPotential = k1Potential(quantile, delta)
		}

		accumulatedSize += nextCentroid.size

	}

	return compressedCentroids

}

func (d *TDigest) quantile(quantile float64) float64 {

	quantileIndex := quantile * float64(d.count())

	maxQuantileIndex := float64(d.centroids[0].size) / 2.0
	if quantileIndex <= maxQuantileIndex {
		return d.centroids[0].mean
	}

	for idx := 0; idx < len(d.centroids)-1; idx++ {

		c1 := d.centroids[idx]
		c2 := d.centroids[idx+1]

		interval := float64(c1.size+c2.size) / 2.0
		if quantileIndex <= maxQuantileIndex+interval {
			k := (quantileIndex - maxQuantileIndex) / interval
			return c1.mean*float64(1-k) + c2.mean*float64(k)
		}
		maxQuantileIndex += interval
	}
	lastCentroid := last(&d.centroids)
	return lastCentroid.mean

}
