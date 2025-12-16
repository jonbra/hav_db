param(
    [string]$ViewerDir = "viewer",
    [string]$VendorDir = "viewer/vendor"
)

Write-Host "Vendoring viewer from $ViewerDir to $VendorDir"

if (-not (Test-Path $ViewerDir)){
    Write-Error "Viewer directory '$ViewerDir' not found"
    exit 2
}

# Ensure vendor dir exists
if (-not (Test-Path $VendorDir)){
    New-Item -ItemType Directory -Path $VendorDir | Out-Null
}

# Copy production build artifacts: .next and public
$srcNext = Join-Path $ViewerDir '.next'
$srcPublic = Join-Path $ViewerDir 'public'

if (Test-Path $srcNext){
    $destNext = Join-Path $VendorDir '.next'
    if (Test-Path $destNext){ Remove-Item -Recurse -Force $destNext }
    Write-Host "Copying $srcNext -> $destNext"
    Copy-Item -Recurse -Force -Path $srcNext -Destination $destNext
} else { Write-Warning ".next build not found in $ViewerDir" }

if (Test-Path $srcPublic){
    $destPublic = Join-Path $VendorDir 'public'
    if (Test-Path $destPublic){ Remove-Item -Recurse -Force $destPublic }
    Write-Host "Copying $srcPublic -> $destPublic"
    Copy-Item -Recurse -Force -Path $srcPublic -Destination $destPublic
} else { Write-Warning "public not found in $ViewerDir" }

Write-Host "Vendoring complete. Review viewer/vendor/ for static assets and server files."
