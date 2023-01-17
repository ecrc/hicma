pipeline {
    agent { label 'jenkinsfile' }
    triggers {
        pollSCM('H/10 * * * *')
    }
    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }
    stages {
        stage('Build_MPI_ON') {
             steps {
                echo "MPI ON :::: BUILD"
                sh '''#!/bin/bash -le

module purge

module load ecrc-extras
module load mkl/2020.0.166
module load gcc/10.2.0
module load cmake/3.19.2
module load hwloc/2.4.0-gcc-10.2.0
module load openmpi/4.1.0-gcc-10.2.0
module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0

module list
# BASH verbose mode
set -x
HICMADEVDIR=$PWD
INSTALLDIR=$PWD/dependencies-prefix
mkdir -p $INSTALLDIR
rm -rf $INSTALLDIR/*
export PKG_CONFIG_PATH=$INSTALLDIR/lib/pkgconfig:$PKG_CONFIG_PATH

# STARS-H
cd $HICMADEVDIR
cd stars-h
git log -1
ls -l
rm -rf build
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DGSL=OFF
make clean
make -j 4
make install

# HCORE
cd $HICMADEVDIR
cd hcore
git log -1
ls -l
rm -rf build
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALLDIR
make clean
make -j 4
make install

# HICMA
cd $HICMADEVDIR
git log -1
ls -l
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=ON
make clean
make -j 4
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
'''

                        stash name: "build-mpion", includes: "build/**"
                    }
                }
                stage('Build_MPI_OFF') {
                    steps {
                        echo "MPI OFF :::: BUILD"
                        sh '''#!/bin/bash -le

module purge

module load ecrc-extras
module load mkl/2020.0.166
module load gcc/10.2.0
module load cmake/3.19.2
module load hwloc/2.4.0-gcc-10.2.0
module load openmpi/4.1.0-gcc-10.2.0
module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0

module list
# BASH verbose mode
set -x
HICMADEVDIR=$PWD
INSTALLDIR=$PWD/dependencies-prefix
mkdir -p $INSTALLDIR
rm -rf $INSTALLDIR/*
export PKG_CONFIG_PATH=$INSTALLDIR/lib/pkgconfig:$PKG_CONFIG_PATH
# STARS-H
cd stars-h
git log -1
ls -l
rm -rf build
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DGSL=OFF
make clean
make -j 4
make install
# HCORE
cd $HICMADEVDIR
cd hcore
git log -1
ls -l
rm -rf build
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALLDIR
make clean
make -j 4
make install
# HICMA
cd $HICMADEVDIR
git log -1
ls -l
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=OFF
make clean
make -j 4
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
'''
                        stash name: "build-mpioff", includes: "build/**"
                    }
                }
                stage('Test_MPI_OFF') {
                    steps {
                        unstash 'build-mpioff'
                        echo "MPI OFF :::: TESTS"
                        sh '''#!/bin/bash -le

module purge

module load ecrc-extras
module load mkl/2020.0.166
module load gcc/10.2.0
module load cmake/3.19.2
module load hwloc/2.4.0-gcc-10.2.0
module load openmpi/4.1.0-gcc-10.2.0
module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0

module list
cd build
rm -rf Testing
ctest -T Test --no-compress-output -V

'''
                    }
                }
        stage('Test_MPI_ON') {
            steps {
                unstash 'build-mpion'
                echo "MPI ON :::: TESTS"
                sh '''#!/bin/bash -le

module purge

module load ecrc-extras
module load mkl/2020.0.166
module load gcc/10.2.0
module load cmake/3.19.2
module load hwloc/2.4.0-gcc-10.2.0
module load openmpi/4.1.0-gcc-10.2.0
module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0

module list
cd build
rm -rf Testing
ctest -T Test --no-compress-output -V

'''
                    }
                }
        stage('Documentation') {
            steps {
                unstash 'build-mpion'
                sh '''#!/bin/bash -el
module purge
module load ecrc-extras
module load mkl/2020.0.166
module load gcc/10.2.0
module load cmake/3.19.2
module load hwloc/2.4.0-gcc-10.2.0
module load openmpi/4.1.0-gcc-10.2.0
module load starpu/1.3.9-gcc-10.2.0-mkl-openmpi-4.1.0

module list
cd $WORKSPACE/build
make docs
'''
                sh '''#!/bin/bash -ex
cd $WORKSPACE
rm -rf   cppcheckhtml
cppcheck --enable=all --xml --xml-version=2 aux/ compute/ control/ hcore/ include/ runtime/ testing/ timing/ -I include/ 2> cppcheck.xml
cppcheck-htmlreport --source-encoding="iso8859-1" --title="HiCMA" --source-dir=. --report-dir=cppcheckhtml --file=cppcheck.xml
'''
                publishHTML( target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false, reportDir: 'build/docs/build/html', reportFiles: 'index.html', reportName: 'Doxygen Documentation', reportTitles: ''] )
                publishHTML( target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false, reportDir: 'cppcheckhtml', reportFiles: 'index.html', reportName: 'CppCheckReport', reportTitles: ''] )
            }
        }
    } 
    // Post build actions
       post {
        //always {
        //}
        //success {
        //}
        //unstable {
        //}
        //failure {
        //}
         unstable {
             emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [culprits(),requestor()]
         }
         failure {
             emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [culprits(),requestor()]
         }
     }
}
