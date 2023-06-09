//////////////////////////////////////
// AMI for AWS Batch (Nextflow)  
// Batch image with 1000G root storage for NF
//////////////////////////////////////

# TODO: install AWS CLI to that AMI


resource "aws_network_interface" "base_batch" {
  subnet_id   = aws_default_subnet.default.id
  tags = {
    Name = "default_network_interface"
  }
}



// Find the latest ECS optimized AMI
data "aws_ami" "ecs_ami" {
  most_recent      = true
  owners           = ["amazon"]

  filter {
    name   = "name"
    values = ["amzn2-ami-ecs-hvm-*"]
  }

  filter {
    name = "architecture"
    values = ["x86_64"]
  }

}

// create a SPOT instance from ami, add 1000 ebs, which is required by Nextflow
resource "aws_instance" "base_batch_nf" {
  ami           = data.aws_ami.ecs_ami.id
  instance_type = "t2.micro"
  //key_name      = "<KEY.PAIR>"

  tags = {
    Name = "spot-base-batch-ami"
  }

  network_interface {
    network_interface_id = aws_network_interface.base_batch.id
    device_index = 0
  }

  root_block_device {
    volume_size = 1000
  }
}

// and now we create ami from the instance above
// the creation of an AMI modelled after an existing EBS-backed EC2 instance
resource "aws_ami_from_instance" "base_batch_nf_ami" {
  name               = "base_batch_nf_ami"
  source_instance_id = aws_instance.base_batch_nf.id
}
