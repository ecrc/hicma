pipeline {
/*
 * Defining where to run
 */
//// Any:
// agent any
//// By agent label:
//      agent { label 'sandybridge' }

    // no agents, each stage must declare it
    agent none
    triggers {
        pollSCM('H/10 * * * *')
    }
    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }


    stages {
        stage('Parallel Build Stage') {
            //when {
            //    branch 'master'
            //}
            // failFast true // abort when one stage fails
            // failFast false // continue even when one stage fails
            failFast false
            parallel {
                stage('MPI ON') {
                    agent {
                        label "jenkinsfile"
                    }
                    steps {
                        echo "MPI ON :::: BUILD"
                        sh '''#!/bin/bash -le

module load mkl/2018-initial
module load gcc/5.5.0
module load cmake/3.9.6
module load hwloc/1.11.8-gcc-5.5.0
module load openmpi/3.0.0-gcc-5.5.0
module load ecrc-extras
module load starpu/1.2.4-gcc-5.5.0-mkl-openmpi-3.0.0
module load gsl/2.4-gcc-5.5.0


module list

# BASH verbose mode
set -x


HICMADEVDIR=$PWD

# Update submodules
#git submodule update --init --recursive

# STARS-H
cd stars-h
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# CHAMELEON
cd $HICMADEVDIR
cd chameleon
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCHAMELEON_USE_MPI=ON  -DCMAKE_INSTALL_PREFIX=$PWD/installdir
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# HICMA
cd $HICMADEVDIR
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=ON
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

'''

                        stash name: "build-mpion", includes: "build/**"
                    }
                }
                stage('MPI OFF') {
                    agent {
                        label "jenkinsfile"
                    }
                    steps {
                        echo "MPI OFF :::: BUILD"
                        sh '''#!/bin/bash -le

module load mkl/2018-initial
module load gcc/5.5.0
module load cmake/3.9.6
module load hwloc/1.11.8-gcc-5.5.0
module load ecrc-extras
module load starpu/1.2.4-gcc-5.5.0-mkl-openmpi-3.0.0
module load gsl/2.4-gcc-5.5.0

module list

# BASH verbose mode
set -x


HICMADEVDIR=$PWD

# Update submodules
#git submodule update --init --recursive

# STARS-H
cd stars-h
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# CHAMELEON
cd $HICMADEVDIR
cd chameleon
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCHAMELEON_USE_MPI=OFF  -DCMAKE_INSTALL_PREFIX=$PWD/installdir
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# HICMA
cd $HICMADEVDIR
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=OFF
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

'''
                        stash name: "build-mpioff", includes: "build/**"
                    }
                }
            }
        }
        stage('Parallel Test Stage') {
            parallel {
                stage('MPI ON') {
                    agent { label 'jenkinsfile' }
                    steps {
                        unstash 'build-mpion'
                        echo "MPI ON :::: TESTS"
                        sh '''#!/bin/bash -le

module load mkl/2018-initial
module load gcc/5.5.0
module load cmake/3.9.6
module load hwloc/1.11.8-gcc-5.5.0
module load openmpi/3.0.0-gcc-5.5.0
module load ecrc-extras
module load starpu/1.2.4-gcc-5.5.0-mkl-openmpi-3.0.0
module load gsl/2.4-gcc-5.5.0

module list

cd build
rm -rf Testing
ctest -T Test --no-compress-output -V
if [ "$BRANCH_NAME" != "master" ]
then
    # valgring for timing takes too long. Test the rest
    #ctest -T memcheck -LE timing 
    :
fi

cd installdir
tar -zcf $WORKSPACE/hicma-mpi-on.tgz ./*
'''
                        step([$class: 'XUnitBuilder',
                             thresholds: [[$class: 'FailedThreshold', unstableThreshold: '0']],
                             tools: [[$class: 'CTestType', pattern: 'build/Testing/**/Test.xml']]])
                        archiveArtifacts allowEmptyArchive: true, artifacts: '*.tgz'
                    }
                }
                stage('MPI OFF') {
                    agent { label 'jenkinsfile' }
                    steps {
                        unstash 'build-mpioff'
                        echo "MPI OFF :::: TESTS"
                        sh '''#!/bin/bash -le

module load mkl/2018-initial
module load gcc/5.5.0
module load cmake/3.9.6
module load hwloc/1.11.8-gcc-5.5.0
module load ecrc-extras
module load starpu/1.2.4-gcc-5.5.0-mkl-openmpi-3.0.0
module load gsl/2.4-gcc-5.5.0

module list

cd build
rm -rf Testing
ctest -T Test --no-compress-output -V
if [ "$BRANCH_NAME" != "master" ]
then
    # valgring for timing takes too long. Test the rest
    #ctest -T memcheck -LE timing 
    :
fi

cd installdir
tar -zcf $WORKSPACE/hicma-mpi-off.tgz ./*
'''
                        step([$class: 'XUnitBuilder',
                             thresholds: [[$class: 'FailedThreshold', unstableThreshold: '0']],
                             tools: [[$class: 'CTestType', pattern: 'build/Testing/**/Test.xml']]])
                        archiveArtifacts allowEmptyArchive: true, artifacts: '*.tgz'
                    }
                }
            }
        }
        stage('Documentation') {
            agent { label 'jenkinsfile' }
            steps {
                unstash 'build-mpion'
                sh '''#!/bin/bash -el
module purge
module load mkl/2018-initial
module load gcc/5.5.0
module load cmake/3.9.6
module load hwloc/1.11.8-gcc-5.5.0
module load openmpi/3.0.0-gcc-5.5.0
module load ecrc-extras
module load starpu/1.2.4-gcc-5.5.0-mkl-openmpi-3.0.0
module load gsl/2.4-gcc-5.5.0
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
             emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
         }
         failure {
             emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
         }
     }
}
