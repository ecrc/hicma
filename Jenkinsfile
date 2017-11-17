pipeline {
/*
 * Defining where to run
 */
//// Any:
// agent any
//// By agent label:
//      agent { label 'sandybridge' }

    agent { label 'Almaha' }
    triggers {
        pollSCM('H/10 * * * *')
    }
    environment {
        XX="gcc"
    }
    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }

    stages {
        stage ('build') {
            steps {
                sh "scripts/build.sh"
            }
        }
        stage ('test') {
            steps {
                sh "echo Spatial Statistics with sqr exp kernel"
                sh "scripts/test.sh 2 - 4 4 --ss"
                sh "echo Electrodynamics with sin kernel"
                sh "scripts/test.sh 1 - 4 4 --edsin"
            }
        }
        stage ('docs') {
            steps {
                sh "cd $WORKSPACE/build && make docs"
                publishHTML( target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false, reportDir: 'build/docs/build/html', reportFiles: 'index.html', reportName: 'Doxygen Documentation', reportTitles: ''] )
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

