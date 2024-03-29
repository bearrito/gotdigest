package gotdigest

import (
	"math"
	"slices"

	"github.com/vmihailenco/msgpack/v5"
)

type Centroid struct {
	Size int
	Mean float64
}

type TDigest struct {
	Centroids     []Centroid
	delta         float64   `msgpack:"-"`
	shouldBuffer  bool      `msgpack:"-"`
	bufferMaxSize int       `msgpack:"-"`
	buffer        []float64 `msgpack:"-"`
	useNormalized bool      `msgpack:"-"`
}

func (d *TDigest) ToMsgPack() ([]byte, error) {

	b, err := msgpack.Marshal(d)
	if err != nil {
		return nil, err

	}
	return b, nil
}

func FromMsgPack(b []byte) (TDigest, error) {

	var digest TDigest
	err := msgpack.Unmarshal(b, &digest)
	if err != nil {
		return digest, err
	}
	return digest, nil
}

func last(c *[]Centroid) *Centroid {

	return &((*c)[len(*c)-1])
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
	digest := TDigest{Centroids: centroids}
	for val := range vals {
		digest.Centroids = append(digest.Centroids, NewCentroidWithValue(float64(val)))

	}

	return digest

}

func (c *Centroid) weight() float64 {

	return c.Mean * float64(c.Size)
}

func (c *Centroid) update(other Centroid) {

	weight := c.weight() + other.weight()
	c.Size += other.Size
	c.Mean = weight / float64(c.Size)
}

func mergeCentroids(first Centroid, second Centroid) Centroid {

	size := (first.Size + second.Size)
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

	return TDigest{Centroids: bins}

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

			digest := TDigest{Centroids: centroids}
			d.merge(&digest)
			d.Centroids = compressBins(d.Centroids, d.delta)
		}

	} else {
		b := NewCentroidWithValue(value)
		digest := NewDigestFromBin(b)
		d.merge(&digest)
		d.Centroids = compressBins(d.Centroids, d.delta)
	}

}

func (d *TDigest) count() int {

	count := 0
	for _, bin := range d.Centroids {
		count += bin.Size

	}
	return count

}

func (d *TDigest) merge(other *TDigest) {

	size := len(d.Centroids) + len(other.Centroids)
	merged := make([]Centroid, 0, size)

	dCount := 0
	otherCount := 0
	for dCount < len(d.Centroids) && otherCount < len(other.Centroids) {

		if d.Centroids[dCount].Mean < other.Centroids[otherCount].Mean {
			merged = append(merged, d.Centroids[dCount])
			dCount++
		} else {

			merged = append(merged, other.Centroids[otherCount])
			otherCount++

		}

	}

	for dCount < len(d.Centroids) {
		merged = append(merged, d.Centroids[dCount])
		dCount++
	}

	for otherCount < len(other.Centroids) {
		merged = append(merged, other.Centroids[otherCount])
		otherCount++
	}

	d.Centroids = merged
}

func compressBins(centroids []Centroid, delta float64) []Centroid {

	if len(centroids) == 0 {
		return centroids
	}

	totalSize := 0
	for _, centroid := range centroids {
		totalSize += centroid.Size

	}

	compressedCentroids := make([]Centroid, 0)
	compressedCentroids = append(compressedCentroids, centroids[0])

	accumulatedSize := compressedCentroids[0].Size

	minPotential := k1Potential(0, delta)
	i := 1
	for i < len(centroids) {
		nextCentroid := centroids[i]
		i++

		quotientIndex := float64((accumulatedSize + nextCentroid.Size)) / float64(totalSize)
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

		accumulatedSize += nextCentroid.Size

	}

	return compressedCentroids

}

func (d *TDigest) quantile(quantile float64) float64 {

	quantileIndex := quantile * float64(d.count())

	maxQuantileIndex := float64(d.Centroids[0].Size) / 2.0
	if quantileIndex <= maxQuantileIndex {
		return d.Centroids[0].Mean
	}

	for idx := 0; idx < len(d.Centroids)-1; idx++ {

		c1 := d.Centroids[idx]
		c2 := d.Centroids[idx+1]

		interval := float64(c1.Size+c2.Size) / 2.0
		if quantileIndex <= maxQuantileIndex+interval {
			k := (quantileIndex - maxQuantileIndex) / interval
			return c1.Mean*float64(1-k) + c2.Mean*float64(k)
		}
		maxQuantileIndex += interval
	}
	lastCentroid := last(&d.Centroids)
	return lastCentroid.Mean

}
