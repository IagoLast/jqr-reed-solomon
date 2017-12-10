const assert = require('assert');
const GF256 = require('jqr-gf');
const ReedSolomon = require('../src/reed-solomon');

describe('reed-solomon', () => {
  const reedSolomon = new ReedSolomon(GF256);
  describe('.decode', () => {
    it('should decode a message with no errors', () => {
      const message = [64, 69, 70, 87, 55, 64, 236, 17, 236, 17, 236, 17, 236, 17, 236, 17, 233, 234, 27, 99, 188, 204, 151, 48, 90, 104];
      const twoS = 10;

      const actual = reedSolomon.decode(message, twoS);

      assert.deepEqual(actual, message);
    });

    it('should fix a message with errors', () => {
      const expected = [64, 69, 70, 87, 55, 64, 236, 17, 236, 17, 236, 17, 236, 17, 236, 17, 233, 234, 27, 99, 188, 204, 151, 48, 90, 104];
      const message = [0, 0, 0, 0, 55, 64, 236, 20, 236, 17, 236, 17, 236, 17, 236, 17, 233, 234, 27, 99, 188, 204, 151, 48, 90, 104];
      const twoS = 10;

      const actual = reedSolomon.decode(message, twoS);

      assert.deepEqual(actual, expected);
    });
  });
});
