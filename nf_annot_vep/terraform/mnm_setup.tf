// terraform init
// terraform plan
// terraform apply

terraform {
  required_version = ">= 0.13"

  // if you need TForm to save its state in S3, uncomment
  // need aws --profile / sh
  /*
  backend "s3" {
    profile = "<AWS_ACCOUNT_NAME>"
    region  = "eu-central-1"
    bucket  = "terraform.<AWS_ACCOUNT_NAME>"
    key     = "aws-env.tfstate"
  }
  */
}

provider "aws" {
  region              = "eu-central-1"
  profile             = "default"
  allowed_account_ids = [<TBD>]
}

