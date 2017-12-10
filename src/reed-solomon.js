const Polynomial = require('jqr-poly').Polynomial;
const GF256 = require('jqr-gf');

class ReedSolomonDecoder {

  constructor() {
    // Assume a QR CODE in a future the field could be passed as parameter
    this.field = GF256; 
  }

  decode(message, twoS) {
    var dataMatrix = false;
    const syndrome = this.computeSyndromes(message, twoS);

    var sigmaOmega = this.runEuclideanAlgorithm(this.field.buildMonomial(twoS, 1), syndrome, twoS);
    var sigma = sigmaOmega[0];
    var omega = sigmaOmega[1];

    var errorLocations = this.findErrorLocations(sigma);
    var errorMagnitudes = this.findErrorMagnitudes(omega, errorLocations, dataMatrix);

    for (let i = 0; i < errorLocations.length; i++) {
      var position = message.length - 1 - this.field.log(errorLocations[i]);
      if (position < 0) {
        throw 'ReedSolomonException Bad error location';
      }
      message[position] = this.field.add(message[position], errorMagnitudes[i]);
    }

    return message;

  }

  computeSyndromes(message, twoS) {
    const syndromeCoefficients = new Array(twoS).fill(0);
    const poly = new Polynomial(message, this.field);
    
    // Calculate syndrome coefficients
    for (let i = 0; i < twoS; i++) {
      var _eval = poly.evaluate(this.field.exp(i));
      syndromeCoefficients[syndromeCoefficients.length - 1 - i] = _eval;
    }
    return new Polynomial(syndromeCoefficients, this.field);
  }

  runEuclideanAlgorithm(a, b, R) {
    // Assume a's degree is >= b's
    if (a.getDegree() < b.getDegree()) {
      var temp = a;
      a = b;
      b = temp;
    }

    var rLast = a;
    var r = b;
    var sLast = this.field.one();
    var s = this.field.zero();
    var tLast = this.field.zero();
    var t = this.field.one();

    // Run Euclidean algorithm until r's degree is less than R/2
    while (r.getDegree() >= Math.floor(R / 2)) {
      var rLastLast = rLast;
      var sLastLast = sLast;
      var tLastLast = tLast;
      rLast = r;
      sLast = s;
      tLast = t;

      // Divide rLastLast by rLast, with quotient in q and remainder in r
      if (rLast.isZero()) {
        // Oops, Euclidean algorithm already terminated?
        throw 'r_{i-1} was zero';
      }
      r = rLastLast;
      var q = this.field.zero();
      var denominatorLeadingTerm = rLast.getCoefficient(rLast.getDegree());
      var dltInverse = this.field.inv(denominatorLeadingTerm);
      while (r.getDegree() >= rLast.getDegree() && !r.isZero()) {
        var degreeDiff = r.getDegree() - rLast.getDegree();
        var scale = this.field.mul(r.getCoefficient(r.getDegree()), dltInverse);
        q = q.add(this.field.buildMonomial(degreeDiff, scale));
        r = r.add(rLast.multiplyByMonomial(degreeDiff, scale));
      }

      s = q.multiplyPolynomial(sLast);
      s = s.add(sLastLast);

      t = q.multiplyPolynomial(tLast);
      t = t.add(tLastLast);
    }

    var sigmaTildeAtZero = t.getCoefficient(0);
    if (sigmaTildeAtZero == 0) {
      throw 'ReedSolomonException sigmaTilde(0) was zero';
    }

    var inverse = this.field.inv(sigmaTildeAtZero);
    var sigma = t.multiplyScalar(inverse);
    var omega = r.multiplyScalar(inverse);
    return [sigma, omega];
  }

  findErrorLocations(errorLocator) {
    // This is a direct application of Chien's search
    var numErrors = errorLocator.getDegree();
    if (numErrors == 1) {
      // shortcut
      return new Array(errorLocator.getCoefficient(1));
    }
    var result = new Array(numErrors);
    var e = 0;
    for (var i = 1; i < 256 && e < numErrors; i++) {
      if (errorLocator.evaluate(i) === 0) {
        result[e] = this.field.inv(i);
        e++;
      }
    }
    if (e != numErrors) {
      throw 'Error locator degree does not match number of roots';
    }
    return result;
  }

  findErrorMagnitudes(errorEvaluator, errorLocations, dataMatrix) {
    // This is directly applying Forney's Formula
    var s = errorLocations.length;
    var result = new Array(s);
    for (var i = 0; i < s; i++) {
      var xiInverse = this.field.inv(errorLocations[i]);
      var denominator = 1;
      for (var j = 0; j < s; j++) {
        if (i != j) {
          denominator = this.field.mul(denominator, GF256.add(1, this.field.mul(errorLocations[j], xiInverse)));
        }
      }
      result[i] = this.field.mul(errorEvaluator.evaluate(xiInverse), this.field.inv(denominator));
      // Thanks to sanfordsquires for this fix:
      if (dataMatrix) {
        result[i] = this.field.mul(result[i], xiInverse);
      }
    }
    return result;
  }

}

module.exports = ReedSolomonDecoder;
