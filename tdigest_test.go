package gotdigest

import (
	"math"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestUtils(t *testing.T) {

	centroids := make([]Centroid, 2)
	centroids[0] = Centroid{1, 10}
	second := Centroid{1, 10}
	centroids[1] = second
	require.Equal(t, second, last(&centroids))

}

func TestK1Potentials(t *testing.T) {

	testCases := []struct {
		title     string
		delta     float64
		potential unNormalizedPotential
	}{
		{
			title:     "Small Delta with K1",
			delta:     5.0,
			potential: k1Potential,
		},
		{
			title:     "Medium Delta with K1",
			delta:     10.0,
			potential: k1Potential,
		},
		{
			title:     "Large Delta with K1",
			delta:     100.0,
			potential: k1Potential,
		},
	}

	for _, tc := range testCases {
		tc := tc
		t.Run(tc.title, func(t *testing.T) {
			delta := tc.delta
			maxPotential := delta / 4.0
			minPotential := -maxPotential

			actualMin := tc.potential(0.0, tc.delta)
			actualMax := tc.potential(1.0, tc.delta)
			require.Equal(t, actualMax, maxPotential)
			require.Equal(t, actualMin, minPotential)

			// Checks Symmetry
			a := tc.potential(0.99, tc.delta)
			b := tc.potential(0.01, tc.delta)
			require.Greater(t, a, b)
			require.Equal(t, a+b, 0.0)
			require.Greater(t, a-b, 1.0)

			require.Equal(t, tc.potential(0.5, tc.delta), 0.0)
			require.Greater(t, b, actualMin)
			require.Greater(t, actualMax, b)

			q8 := tc.potential(0.8, tc.delta)
			q2 := tc.potential(0.2, tc.delta)
			require.GreaterOrEqual(t, q8-q2, 1.0)

		})
	}
}

func TestK0Potentials(t *testing.T) {

	testCases := []struct {
		title     string
		delta     float64
		potential unNormalizedPotential
	}{
		{
			title:     "Small Delta with K0",
			delta:     5.0,
			potential: k0Potential,
		},
		{
			title:     "Medium Delta with K0",
			delta:     10.0,
			potential: k0Potential,
		},
		{
			title:     "Large Delta with K0",
			delta:     100.0,
			potential: k0Potential,
		},
	}

	for _, tc := range testCases {
		tc := tc
		t.Run(tc.title, func(t *testing.T) {
			delta := tc.delta
			maxPotential := delta / 2.0
			minPotential := 0.0

			actualMin := tc.potential(0.0, tc.delta)
			actualMax := tc.potential(1.0, tc.delta)
			require.Equal(t, actualMax, maxPotential)
			require.Equal(t, actualMin, minPotential)

			// Checks Symmetry
			a := tc.potential(0.99, tc.delta)
			b := tc.potential(0.01, tc.delta)
			require.Greater(t, a, b)
			require.Greater(t, actualMax, a)
			require.Greater(t, b, actualMin)

			require.Equal(t, tc.potential(0.5, tc.delta), (delta/2.0)*0.5)

			q8 := tc.potential(0.8, tc.delta)
			q2 := tc.potential(0.2, tc.delta)
			require.GreaterOrEqual(t, q8-q2, 1.0)

		})
	}
}

func TestK2Potentials(t *testing.T) {

	testCases := []struct {
		title     string
		delta     float64
		potential normalizedPotential
		n         int
	}{
		{
			title:     "Small Delta with K0",
			delta:     20.0,
			potential: k2Potential,
			n:         10000,
		},
		{
			title:     "Medium Delta with K0",
			delta:     100.0,
			potential: k2Potential,
			n:         10000,
		},
		{
			title:     "Large Delta with K0",
			delta:     500.0,
			potential: k2Potential,
			n:         10000,
		},
	}

	for _, tc := range testCases {
		tc := tc
		t.Run(tc.title, func(t *testing.T) {

			// maxPotential := delta / 2.0
			// minPotential := 0.0

			actualMin := tc.potential(0.0, tc.delta, tc.n)
			actualMax := tc.potential(1.0, tc.delta, tc.n)
			// require.Equal(t, actualMax, maxPotential)
			// require.Equal(t, actualMin, minPotential)

			// Checks Symmetry
			a := tc.potential(0.99, tc.delta, tc.n)
			b := tc.potential(0.01, tc.delta, tc.n)
			require.Greater(t, a, b)
			require.Less(t, a+b, 0.000001)
			require.Greater(t, actualMax, a)
			require.Greater(t, b, actualMin)

			require.Equal(t, 0, 0, tc.potential(0.5, tc.delta, tc.n))

			q8 := tc.potential(0.8, tc.delta, tc.n)
			q2 := tc.potential(0.2, tc.delta, tc.n)
			require.GreaterOrEqual(t, q8-q2, 1.0)

		})
	}
}

func TestBin(t *testing.T) {

	bin := Centroid{1, 10}
	require.Equal(t, bin.Size, 1)
	require.Equal(t, bin.Mean, float64(10))

	require.Equal(t, bin.weight(), float64(10))

	anotherBin := Centroid{2, 10}
	require.Equal(t, anotherBin.weight(), float64(20.0))

	mergedBin := mergeCentroids(bin, anotherBin)

	require.Equal(t, mergedBin.Size, 3)
	require.Equal(t, mergedBin.Mean, float64(10.0))

	bin.update(anotherBin)
	require.Equal(t, bin.Size, 3)
	require.Equal(t, bin.Mean, float64(10.0))

}

func TestTDigestConstructionAndFilling(t *testing.T) {

	digest := NewDigest(10, false, 0, false)
	require.NotNil(t, digest)
	require.Equal(t, digest.count(), 0)

	digest.append(1.0)
	require.Equal(t, digest.count(), 1)
	digest.append(2.0)
	require.Equal(t, digest.count(), 2)
	i := 0
	for i < 100 {

		digest.append(float64(i))
		i++

	}
	require.Equal(t, digest.count(), 102)

	first := float64(1.5)
	second := float64(2.5)
	third := float64(3.5)
	digest = NewDigestFromValues(first, second, third)
	require.NotNil(t, digest)
	require.Equal(t, digest.Centroids[0].Mean, first)
	require.Equal(t, digest.Centroids[1].Mean, second)
	require.Equal(t, digest.Centroids[2].Mean, third)
}

func TestTDigestMerge(t *testing.T) {

	digest1 := NewDigest(10.0, false, 0, false)
	digest2 := NewDigest(10.0, false, 0, false)
	for i := 0; i < 100; i++ {

		// Sorted Set is a TDigest
		digest1.Centroids = append(digest1.Centroids, Centroid{1, float64(i * 2)})
		digest2.Centroids = append(digest2.Centroids, Centroid{1, float64(i*2 + 1)})

	}

	// This only merges but does not compress
	// The result should just be a sorted set
	digest1.merge(&digest2)

	require.Equal(t, len(digest1.Centroids), 200)
	for i := 0; i < 200; i++ {

		require.Equal(t, digest1.Centroids[i].Mean, float64(i))
	}

}

func TestTDigestStatisticsUnNormalized(t *testing.T) {

	testCases := []struct {
		title        string
		delta        float64
		shouldBuffer bool
		bufferMax    int
	}{
		{
			title:        "Small Delta with buffering",
			delta:        10.0,
			shouldBuffer: true,
			bufferMax:    10,
		},
		{
			title: "Small Delta",
			delta: 10.0,
		},
		{
			title: "Medium Delta",
			delta: 100.0,
		},
		{
			title: "Large Delta",
			delta: 500.0},
	}

	for _, tc := range testCases {
		tc := tc
		t.Run(tc.title, func(t *testing.T) {

			delta := tc.delta
			digest := NewDigest(delta, tc.shouldBuffer, tc.bufferMax, false)
			maxData := 0.0
			minData := 0.0
			i := 0
			for i < 10000 {
				data := rand.Float64()*10000 + rand.Float64()*float64(i)
				if i == 0 {
					maxData = data
					minData = data

				} else {
					maxData = math.Max(maxData, data)
					minData = math.Min(minData, data)

				}

				digest.append(data)
				i++
			}

			// require.Greater(t, float64(len(digest.centroids)), math.Floor(delta/2.0))
			require.Less(t, float64(len(digest.Centroids)), math.Ceil(delta))

			minCentroid := digest.quantile(0.0)
			maxCentroid := digest.quantile(1.0)

			require.Greater(t, maxCentroid, minCentroid)
			require.GreaterOrEqual(t, minCentroid, minData)
			require.LessOrEqual(t, maxCentroid, maxData)

			require.Greater(t, digest.quantile(0.1), minData)
			require.Less(t, digest.quantile(0.99), maxData)

			require.Less(t, digest.quantile(0.05), maxCentroid)
			require.Greater(t, digest.quantile(0.05), minCentroid)

			for i := 0; i < len(digest.Centroids)-1; i++ {

				require.Less(t, digest.Centroids[i].Mean, digest.Centroids[i+1].Mean)
			}
		})
	}

}

func TestTDigestStatisticsNormalized(t *testing.T) {

	testCases := []struct {
		title        string
		delta        float64
		shouldBuffer bool
		bufferMax    int
	}{
		{
			title:        "Small Delta with buffering",
			delta:        10.0,
			shouldBuffer: true,
			bufferMax:    10,
		},
		{
			title: "Small Delta",
			delta: 50.0,
		},
		{
			title: "Medium Delta",
			delta: 100.0,
		},
		{
			title: "Large Delta",
			delta: 500.0},
	}

	for _, tc := range testCases {
		tc := tc
		t.Run(tc.title, func(t *testing.T) {

			delta := tc.delta
			digest := NewDigest(delta, tc.shouldBuffer, tc.bufferMax, true)
			maxData := 0.0
			minData := 0.0
			i := 0
			for i < 10000 {
				data := rand.Float64()*10000 + rand.Float64()*float64(i)
				if i == 0 {
					maxData = data
					minData = data

				} else {
					maxData = math.Max(maxData, data)
					minData = math.Min(minData, data)

				}

				digest.append(data)
				i++
			}

			// require.Greater(t, float64(len(digest.centroids)), math.Floor(delta/2.0))
			require.Less(t, float64(len(digest.Centroids)), math.Ceil(delta))

			minCentroid := digest.quantile(0.0)
			maxCentroid := digest.quantile(1.0)

			require.Greater(t, maxCentroid, minCentroid)
			require.GreaterOrEqual(t, minCentroid, minData)
			require.LessOrEqual(t, maxCentroid, maxData)

			require.Greater(t, digest.quantile(0.1), minData)
			require.Less(t, digest.quantile(0.99), maxData)

			require.Less(t, digest.quantile(0.05), maxCentroid)
			require.Greater(t, digest.quantile(0.05), minCentroid)

			for i := 0; i < len(digest.Centroids)-1; i++ {

				require.Less(t, digest.Centroids[i].Mean, digest.Centroids[i+1].Mean)
			}
		})
	}

}

func TestMarshalling(t *testing.T) {

	digest := TDigest{}
	i := 0
	centroids := make([]Centroid, 0, 100)
	for i < 100 {

		centroids = append(centroids, Centroid{1, float64(i)})
		i++
	}
	digest.Centroids = centroids

	b, err := digest.ToMsgPack()
	require.Nil(t, err)
	require.Greater(t, len(b), 0)
	actualDigest, err := FromMsgPack(b)
	require.Nil(t, err)
	require.Equal(t, digest, actualDigest)

}
