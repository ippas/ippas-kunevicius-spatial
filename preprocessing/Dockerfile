FROM rocker/tidyverse:4.5.2

LABEL maintainer="Mateusz Zięba"
LABEL description="RStudio Server and SSH environment for ippas-kunevicius-spatial"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Europe/Warsaw \
    DEFAULT_USER=mateusz

# Replace the default RStudio user with mateusz.
RUN if getent passwd 1000 > /dev/null; then \
        userdel --remove "$(getent passwd 1000 | cut -d: -f1)"; \
    fi \
    && /rocker_scripts/default_user.sh

# Tools needed for R development, package installation and SSH access.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    git \
    less \
    nano \
    sudo \
    openssh-server \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libhdf5-dev \
    libglpk-dev \
    libgit2-dev \
    libicu-dev \
    libsqlite3-dev \
    libxt-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# SSH configuration: key-based login for Positron.
RUN mkdir -p \
        /run/sshd \
        /home/mateusz/.ssh \
        /home/mateusz/projects/ippas-kunevicius-spatial \
    && chown -R mateusz:mateusz /home/mateusz \
    && chmod 700 /home/mateusz/.ssh \
    && printf '%s\n' \
        'PermitRootLogin no' \
        'PasswordAuthentication no' \
        'KbdInteractiveAuthentication no' \
        'PubkeyAuthentication yes' \
        'AuthorizedKeysFile .ssh/authorized_keys' \
        'AllowUsers mateusz' \
        'AllowTcpForwarding yes' \
        > /etc/ssh/sshd_config.d/99-positron.conf

WORKDIR /home/mateusz/projects/ippas-kunevicius-spatial

EXPOSE 8787
EXPOSE 22

# Start SSH in the background, then start RStudio Server.
CMD ["/bin/bash", "-c", "ssh-keygen -A && /usr/sbin/sshd && exec /init"]
