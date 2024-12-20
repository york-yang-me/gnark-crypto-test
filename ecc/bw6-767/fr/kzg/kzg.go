// Copyright 2020 ConsenSys Software Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Code generated by consensys/gnark-crypto DO NOT EDIT

package kzg

import (
	"errors"
	"math/big"
	"runtime"
	"sync"

	"crypto/sha256"
	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bw6-767"
	"github.com/consensys/gnark-crypto/ecc/bw6-767/fr"
	"github.com/consensys/gnark-crypto/ecc/bw6-767/fr/fft"
	"github.com/consensys/gnark-crypto/ecc/bw6-767/fr/polynomial"
	"github.com/consensys/gnark-crypto/fiat-shamir"
)

var (
	ErrInvalidNbDigests              = errors.New("number of digests is not the same as the number of polynomials")
	ErrInvalidPolynomialSize         = errors.New("invalid polynomial size (larger than SRS or == 0)")
	ErrVerifyOpeningProof            = errors.New("can't verify opening proof")
	ErrVerifyBatchOpeningSinglePoint = errors.New("can't verify batch opening proof at single point")
	ErrInvalidDomain                 = errors.New("domain cardinality is smaller than polynomial degree")
	ErrMinSRSSize                    = errors.New("minimum srs size is 2")
)

// Digest commitment of a polynomial.
type Digest = bw6767.G1Affine

// SRS stores the result of the MPC
type SRS struct {
	G1 []bw6767.G1Affine  // [gen [alpha]gen , [alpha**2]gen, ... ]
	G2 [2]bw6767.G2Affine // [gen, [alpha]gen ]
}

// NewSRS returns a new SRS using alpha as randomness source
//
// In production, a SRS generated through MPC should be used.
//
// implements io.ReaderFrom and io.WriterTo
func NewSRS(size uint64, bAlpha *big.Int) (*SRS, error) {
	if size < 2 {
		return nil, ErrMinSRSSize
	}
	var srs SRS
	srs.G1 = make([]bw6767.G1Affine, size)

	var alpha fr.Element
	alpha.SetBigInt(bAlpha)

	_, _, gen1Aff, gen2Aff := bw6767.Generators()
	srs.G1[0] = gen1Aff
	srs.G2[0] = gen2Aff
	srs.G2[1].ScalarMultiplication(&gen2Aff, bAlpha)

	alphas := make([]fr.Element, size-1)
	alphas[0] = alpha
	for i := 1; i < len(alphas); i++ {
		alphas[i].Mul(&alphas[i-1], &alpha)
	}
	for i := 0; i < len(alphas); i++ {
		alphas[i].FromMont()
	}
	g1s := bw6767.BatchScalarMultiplicationG1(&gen1Aff, alphas)
	copy(srs.G1[1:], g1s)

	return &srs, nil
}

// OpeningProof KZG proof for opening at a single point.
//
// implements io.ReaderFrom and io.WriterTo
type OpeningProof struct {
	// H quotient polynomial (f - f(z))/(x-z)
	H bw6767.G1Affine

	// Point at which the polynomial is evaluated
	Point fr.Element

	// ClaimedValue purported value
	ClaimedValue fr.Element
}

// BatchOpeningProof opening proof for many polynomials at the same point
//
// implements io.ReaderFrom and io.WriterTo
type BatchOpeningProof struct {
	// H quotient polynomial Sum_i gamma**i*(f - f(z))/(x-z)
	H bw6767.G1Affine

	// Point at which the polynomials are evaluated
	Point fr.Element

	// ClaimedValues purported values
	ClaimedValues []fr.Element
}

// Commit commits to a polynomial using a multi exponentiation with the SRS.
// It is assumed that the polynomial is in canonical form, in Montgomery form.
//
// if a ecc.CPUSemaphore is NOT provided, this function may split the multi exponentiation
// to use all available CPUs
func Commit(p polynomial.Polynomial, srs *SRS, opt ...*ecc.CPUSemaphore) (Digest, error) {
	if len(p) == 0 || len(p) > len(srs.G1) {
		return Digest{}, ErrInvalidPolynomialSize
	}

	var res bw6767.G1Affine

	config := ecc.MultiExpConfig{ScalarsMont: true}
	if len(opt) > 0 {
		config.CPUSemaphore = opt[0]
		if _, err := res.MultiExp(srs.G1[:len(p)], p, config); err != nil {
			return Digest{}, err
		}
	} else {
		numCpus := runtime.NumCPU()
		if numCpus > 16 {
			// we split the multiExp in 2.
			// TODO with machines with more than 32 physical cores, we may want to split in more chunks
			m := len(p) / 2
			chPart1 := make(chan struct{}, 1)
			config.CPUSemaphore = ecc.NewCPUSemaphore(numCpus)
			var p1, p2 bw6767.G1Jac
			var err1 error
			go func() {
				_, err1 = p1.MultiExp(srs.G1[:m], p[:m], config)
				close(chPart1)
			}()
			// part 2
			if _, err := p2.MultiExp(srs.G1[m:len(p)], p[m:], config); err != nil {
				return Digest{}, err
			}
			<-chPart1 // wait for part 1
			if err1 != nil {
				return Digest{}, err1
			}

			p1.AddAssign(&p2)
			res.FromJacobian(&p1)
		} else {
			if _, err := res.MultiExp(srs.G1[:len(p)], p, config); err != nil {
				return Digest{}, err
			}
		}
	}

	return res, nil
}

// Open computes an opening proof of polynomial p at given point.
// fft.Domain Cardinality must be larger than p.Degree()
func Open(p polynomial.Polynomial, point *fr.Element, domain *fft.Domain, srs *SRS) (OpeningProof, error) {
	if len(p) == 0 || len(p) > len(srs.G1) {
		return OpeningProof{}, ErrInvalidPolynomialSize
	}
	if len(p) > int(domain.Cardinality) {
		return OpeningProof{}, ErrInvalidDomain
	}

	// build the proof
	res := OpeningProof{
		Point:        *point,
		ClaimedValue: p.Eval(point),
	}

	// compute H
	_p := make(polynomial.Polynomial, len(p), domain.Cardinality)
	copy(_p, p)
	h := dividePolyByXminusA(domain, _p, res.ClaimedValue, res.Point)
	_p = nil // h re-use this memory

	// commit to H
	hCommit, err := Commit(h, srs)
	if err != nil {
		return OpeningProof{}, err
	}
	res.H.Set(&hCommit)

	return res, nil
}

// Verify verifies a KZG opening proof at a single point
func Verify(commitment *Digest, proof *OpeningProof, srs *SRS) error {

	// comm(f(a))
	var claimedValueG1Aff bw6767.G1Affine
	var claimedValueBigInt big.Int
	proof.ClaimedValue.ToBigIntRegular(&claimedValueBigInt)
	claimedValueG1Aff.ScalarMultiplication(&srs.G1[0], &claimedValueBigInt)

	// [f(alpha) - f(a)]G1Jac
	var fminusfaG1Jac, tmpG1Jac bw6767.G1Jac
	fminusfaG1Jac.FromAffine(commitment)
	tmpG1Jac.FromAffine(&claimedValueG1Aff)
	fminusfaG1Jac.SubAssign(&tmpG1Jac)

	// [-H(alpha)]G1Aff
	var negH bw6767.G1Affine
	negH.Neg(&proof.H)

	// [alpha-a]G2Jac
	var alphaMinusaG2Jac, genG2Jac, alphaG2Jac bw6767.G2Jac
	var pointBigInt big.Int
	proof.Point.ToBigIntRegular(&pointBigInt)
	genG2Jac.FromAffine(&srs.G2[0])
	alphaG2Jac.FromAffine(&srs.G2[1])
	alphaMinusaG2Jac.ScalarMultiplication(&genG2Jac, &pointBigInt).
		Neg(&alphaMinusaG2Jac).
		AddAssign(&alphaG2Jac)

	// [alpha-a]G2Aff
	var xminusaG2Aff bw6767.G2Affine
	xminusaG2Aff.FromJacobian(&alphaMinusaG2Jac)

	// [f(alpha) - f(a)]G1Aff
	var fminusfaG1Aff bw6767.G1Affine
	fminusfaG1Aff.FromJacobian(&fminusfaG1Jac)

	// e([-H(alpha)]G1Aff, G2gen).e([-H(alpha)]G1Aff, [alpha-a]G2Aff) ==? 1
	check, err := bw6767.PairingCheck(
		[]bw6767.G1Affine{fminusfaG1Aff, negH},
		[]bw6767.G2Affine{srs.G2[0], xminusaG2Aff},
	)
	if err != nil {
		return err
	}
	if !check {
		return ErrVerifyOpeningProof
	}
	return nil
}

// BatchOpenSinglePoint creates a batch opening proof at _val of a list of polynomials.
// It's an interactive protocol, made non interactive using Fiat Shamir.
// point is the point at which the polynomials are opened.
// digests is the list of committed polynomials to open, need to derive the challenge using Fiat Shamir.
// polynomials is the list of polynomials to open.
func BatchOpenSinglePoint(polynomials []polynomial.Polynomial, digests []Digest, point *fr.Element, domain *fft.Domain, srs *SRS) (BatchOpeningProof, error) {

	// check for invalid sizes
	nbDigests := len(digests)
	if nbDigests != len(polynomials) {
		return BatchOpeningProof{}, ErrInvalidNbDigests
	}
	largestPoly := -1
	for _, p := range polynomials {
		if len(p) == 0 || len(p) > len(srs.G1) {
			return BatchOpeningProof{}, ErrInvalidPolynomialSize
		}
		if len(p) > int(domain.Cardinality) {
			return BatchOpeningProof{}, ErrInvalidDomain
		}
		if len(p) > largestPoly {
			largestPoly = len(p)
		}
	}

	var res BatchOpeningProof

	// compute the purported values
	res.ClaimedValues = make([]fr.Element, len(polynomials))
	var wg sync.WaitGroup
	wg.Add(len(polynomials))
	for i := 0; i < len(polynomials); i++ {
		go func(at int) {
			res.ClaimedValues[at] = polynomials[at].Eval(point)
			wg.Done()
		}(i)
	}

	// set the point at which the evaluation is done
	res.Point = *point

	// derive the challenge gamma, binded to the point and the commitments
	fs := fiatshamir.NewTranscript(sha256.New(), "gamma")
	if err := fs.Bind("gamma", res.Point.Marshal()); err != nil {
		return BatchOpeningProof{}, err
	}
	for i := 0; i < len(digests); i++ {
		if err := fs.Bind("gamma", digests[i].Marshal()); err != nil {
			return BatchOpeningProof{}, err
		}
	}
	gammaByte, err := fs.ComputeChallenge("gamma")
	if err != nil {
		return BatchOpeningProof{}, err
	}
	var gamma fr.Element
	gamma.SetBytes(gammaByte)

	// compute sum_i gamma**i*f(a)
	var sumGammaiTimesEval fr.Element
	chSumGammai := make(chan struct{}, 1)
	go func() {
		// wait for polynomial evaluations to be completed (res.ClaimedValues)
		wg.Wait()
		sumGammaiTimesEval = res.ClaimedValues[nbDigests-1]
		for i := nbDigests - 2; i >= 0; i-- {
			sumGammaiTimesEval.Mul(&sumGammaiTimesEval, &gamma).
				Add(&sumGammaiTimesEval, &res.ClaimedValues[i])
		}
		close(chSumGammai)
	}()

	// compute sum_i gamma**i*f
	// that is p0 + gamma * p1 + gamma^2 * p2 + ... gamma^n * pn
	// note: if we are willing to paralellize that, we could clone the poly and scale them by
	// gamma n in parallel, before reducing into sumGammaiTimesPol
	sumGammaiTimesPol := make(polynomial.Polynomial, largestPoly, domain.Cardinality)
	copy(sumGammaiTimesPol, polynomials[0])
	gammaN := gamma
	var pj fr.Element
	for i := 1; i < len(polynomials); i++ {
		for j := 0; j < len(polynomials[i]); j++ {
			pj.Mul(&polynomials[i][j], &gammaN)
			sumGammaiTimesPol[j].Add(&sumGammaiTimesPol[j], &pj)
		}
		// TODO @thomas check that we don't need to scale by gamma if polynmials[i] < largestPoly
		gammaN.Mul(&gammaN, &gamma)
	}

	// compute H
	<-chSumGammai
	h := dividePolyByXminusA(domain, sumGammaiTimesPol, sumGammaiTimesEval, res.Point)
	sumGammaiTimesPol = nil // same memory as h

	res.H, err = Commit(h, srs)
	if err != nil {
		return BatchOpeningProof{}, err
	}

	return res, nil
}

// BatchVerifySinglePoint verifies a batched opening proof at a single point of a list of polynomials.
// point: point at which the polynomials are evaluated
// claimedValues: claimed values of the polynomials at _val
// commitments: list of commitments to the polynomials which are opened
// batchOpeningProof: the batched opening proof at a single point of the polynomials.
func BatchVerifySinglePoint(digests []Digest, batchOpeningProof *BatchOpeningProof, srs *SRS) error {
	nbDigests := len(digests)

	// check consistancy between numbers of claims vs number of digests
	if nbDigests != len(batchOpeningProof.ClaimedValues) {
		return ErrInvalidNbDigests
	}

	// derive the challenge gamma, binded to the point and the commitments
	fs := fiatshamir.NewTranscript(sha256.New(), "gamma")
	if err := fs.Bind("gamma", batchOpeningProof.Point.Marshal()); err != nil {
		return err
	}
	for i := 0; i < len(digests); i++ {
		if err := fs.Bind("gamma", digests[i].Marshal()); err != nil {
			return err
		}
	}
	gammaByte, err := fs.ComputeChallenge("gamma")
	if err != nil {
		return err
	}
	var gamma fr.Element
	gamma.SetBytes(gammaByte)

	var sumGammaiTimesEval fr.Element
	sumGammaiTimesEval.Set(&batchOpeningProof.ClaimedValues[nbDigests-1])
	for i := nbDigests - 2; i >= 0; i-- {
		sumGammaiTimesEval.Mul(&sumGammaiTimesEval, &gamma).
			Add(&sumGammaiTimesEval, &batchOpeningProof.ClaimedValues[i])
	}

	var sumGammaiTimesEvalBigInt big.Int
	sumGammaiTimesEval.ToBigIntRegular(&sumGammaiTimesEvalBigInt)
	var sumGammaiTimesEvalG1Aff bw6767.G1Affine
	sumGammaiTimesEvalG1Aff.ScalarMultiplication(&srs.G1[0], &sumGammaiTimesEvalBigInt)

	var acc fr.Element
	acc.SetOne()
	gammai := make([]fr.Element, len(digests))
	gammai[0][0] = 1 // 1 in regular form
	for i := 1; i < len(digests); i++ {
		acc.Mul(&acc, &gamma)
		gammai[i] = acc
		gammai[i].FromMont()
	}
	var sumGammaiTimesDigestsG1Aff bw6767.G1Affine

	if _, err := sumGammaiTimesDigestsG1Aff.MultiExp(digests, gammai, ecc.MultiExpConfig{}); err != nil {
		return err
	}

	// sum_i [gamma**i * (f-f(a))]G1
	var sumGammiDiffG1Aff bw6767.G1Affine
	var t1, t2 bw6767.G1Jac
	t1.FromAffine(&sumGammaiTimesDigestsG1Aff)
	t2.FromAffine(&sumGammaiTimesEvalG1Aff)
	t1.SubAssign(&t2)
	sumGammiDiffG1Aff.FromJacobian(&t1)

	// [alpha-a]G2Jac
	var alphaMinusaG2Jac, genG2Jac, alphaG2Jac bw6767.G2Jac
	var pointBigInt big.Int
	batchOpeningProof.Point.ToBigIntRegular(&pointBigInt)
	genG2Jac.FromAffine(&srs.G2[0])
	alphaG2Jac.FromAffine(&srs.G2[1])
	alphaMinusaG2Jac.ScalarMultiplication(&genG2Jac, &pointBigInt).
		Neg(&alphaMinusaG2Jac).
		AddAssign(&alphaG2Jac)

	// [alpha-a]G2Aff
	var xminusaG2Aff bw6767.G2Affine
	xminusaG2Aff.FromJacobian(&alphaMinusaG2Jac)

	// [-H(alpha)]G1Aff
	var negH bw6767.G1Affine
	negH.Neg(&batchOpeningProof.H)

	// check the pairing equation
	check, err := bw6767.PairingCheck(
		[]bw6767.G1Affine{sumGammiDiffG1Aff, negH},
		[]bw6767.G2Affine{srs.G2[0], xminusaG2Aff},
	)
	if err != nil {
		return err
	}
	if !check {
		return ErrVerifyBatchOpeningSinglePoint
	}
	return nil
}

// dividePolyByXminusA computes (f-f(a))/(x-a), in canonical basis, in regular form
// f memory is re-used for the result
// @precondition: cap(f) must be d.Cardinality or this will panic
func dividePolyByXminusA(d *fft.Domain, f polynomial.Polynomial, fa, a fr.Element) polynomial.Polynomial {
	degree := len(f) - 1 // the result if of degree  deg(f)-1

	// first we compute f-f(a)
	f = f[:d.Cardinality]
	f[0].Sub(&f[0], &fa)

	// now we use syntetic division to divide by x-a
	// TODO check with large polynomials, we may want to parallelize this
	var c, t fr.Element
	for i := len(f) - 1; i >= 0; i-- {
		t.Mul(&c, &a)
		f[i].Add(&f[i], &t)
		c, f[i] = f[i], c
	}

	// the result is of degree deg(f)-1
	return f[:degree]
}
