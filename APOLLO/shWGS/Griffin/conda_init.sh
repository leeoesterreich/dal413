# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ihome/crc/install/python/ondemand-jupyter-python3.10-2023.03/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ihome/crc/install/python/ondemand-jupyter-python3.10-2023.03/etc/profile.d/conda.sh" ]; then
        . "/ihome/crc/install/python/ondemand-jupyter-python3.10-2023.03/etc/profile.d/conda.sh"
    else
        export PATH="/ihome/crc/install/python/ondemand-jupyter-python3.10-2023.03/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<< 