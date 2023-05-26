///////////////////////////////////////////////////
// entry point instance to run analyses pipelines
// NOTE: need to ss and install conda / terraform over there
///////////////////////////////////////////////////


// find the most recent ubuntu ami information
// https://registry.terraform.io/providers/hashicorp/aws/latest/docs/resources/instance
data "aws_ami" "ubuntu" {
  most_recent = true

  filter {
    name   = "name"
    values = ["ubuntu/images/hvm-ssd/ubuntu-focal-20.04-amd64-server-*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }

  owners = ["099720109477"] # from Canonical owner
}

// launch AMI 
resource "aws_instance" "basic" {
  ami = data.aws_ami.ubuntu.id 
  instance_type = "t2.micro"

  key_name      = "<KEY.PAIR>"  # for ssh user who will connect (create separatelly)
  security_groups = [aws_security_group.allow_ssh.name] # allow traffic for ssh

  root_block_device {
    volume_size = 30 
  }

  ebs_block_device {
    device_name = "/dev/xvdb"
    volume_size = 100
    volume_type = "standard"
  }

  tags = {
    Name = "basic"
  }
}

// use default VPC
resource "aws_default_vpc" "default" {
  tags = {
    Name = "Default VPC"
  }
}

resource "aws_security_group" "allow_ssh" {
  name        = "allow_ssh"
  description = "Allow SSH inbound traffic"
  vpc_id      = aws_default_vpc.default.id 

  ingress {
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = "allow_ssh"
  }
}

