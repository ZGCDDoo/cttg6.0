
make -j 4   FourierTests \
            FillingAndDoccTests \
            GreenBinningTests \
            GreenMatTests \
            GreenTauTests \
            H0TriangleTests \
            IntegratorTests \
            ModelTriangle2x2Tests \
            UtilitiesTests \
            SelfConsistencyTests \
            MarkovChainTests \
            MarkovChainSquare2x2Tests \
            MarkovChainAuxTests \
            MarkovChainSubMatrixSquare2x2Tests \
            MarkovChainAuxSubMatrixTests \
            MarkovChainSubMatrixTests \
            MatrixTests \
            ObservablesTests ;

make test
