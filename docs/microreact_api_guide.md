# Microreact API Integration Guide

This guide explains how to use the Microreact API integration in the HAV Database application to upload phylogenetic analyses directly to [microreact.org](https://microreact.org).

## Table of Contents

- [Overview](#overview)
- [Getting Started](#getting-started)
- [Uploading Datasets](#uploading-datasets)
- [Team Management](#team-management)
- [Viewing Projects](#viewing-projects)
- [Troubleshooting](#troubleshooting)
- [API Reference](#api-reference)

---

## Overview

The Microreact integration allows you to:

- **Upload analyses directly** to microreact.org from the application
- **Upload existing .microreact files** without re-running analyses
- **Create and manage teams** for collaborative viewing
- **Share projects** with team members automatically after upload
- **Track uploaded projects** within the session

## Getting Started

### Step 1: Create a Microreact Account

1. Go to [microreact.org](https://microreact.org)
2. Click **Sign In** and create an account (you can use Google, GitHub, or email)
3. Verify your email if required

### Step 2: Obtain Your API Access Token

1. Log in to [microreact.org](https://microreact.org)
2. Navigate to **My Account** → **Account Settings**  
   Direct link: [https://microreact.org/my-account/settings](https://microreact.org/my-account/settings)
3. Find the **API Access** section
4. Copy your **Access Token** (it starts with `eyJ...`)

> ⚠️ **Important**: Keep your API token secure. Anyone with your token can create projects under your account.

### Step 3: Configure the Token in the Application

1. Open the HAV Database application
2. Go to the **Settings** tab
3. In the **Microreact API Settings** box:
   - Paste your API token into the **API Token** field
   - Click **Save Token**
   - Optionally click **Test Connection** to verify the format

Your token is now saved for this session. You'll need to re-enter it when you restart the application.

---

## Uploading Datasets

### Option 1: Upload Current Analysis

This uploads your current phylogenetic analysis (tree + metadata) to Microreact.

**Prerequisites:**
- Complete an alignment and build a tree in the **Analysis** tab
- Have your API token configured in **Settings**

**Steps:**

1. Go to the **Microreact** tab
2. Enter a **Project Name** (default includes today's date)
3. Enter a **Description** for your project
4. Click **Upload to Microreact Server**
5. Wait for the upload to complete
6. The project URL will appear in the "Microreact Visualization" panel

### Option 2: Upload an Existing .microreact File

If you have a previously saved `.microreact` file, you can upload it directly.

**Steps:**

1. Go to the **Microreact** tab
2. In the **Upload .microreact File** box:
   - Click **Choose .microreact file** and select your file
   - Click **Upload File to Server**
3. The project URL will appear once uploaded

### Option 3: Manual Upload (No API Required)

If you prefer not to use the API, you can download files and upload manually:

1. Click **Download Files (ZIP)** to get `metadata.csv` and `tree.nwk`
2. Go to [microreact.org/upload](https://microreact.org/upload)
3. Drag and drop both files
4. Configure your visualization and save

---

## Team Management

Teams allow you to share projects with a group of collaborators. This is an **experimental feature** from Microreact.

### Creating a Team

1. Go to the **Settings** tab
2. In the **Team Management** section, find **Create New Team**
3. Enter a **Team Name**
4. Click **Create Team**
5. The Team ID will be saved automatically

### Using an Existing Team

If you already have a team ID:

1. Go to the **Settings** tab
2. Enter the **Team ID** in the text field
3. Click **Save Team ID**

### Managing Team Members

#### Adding Members

1. In **Settings** → **Team Members**
2. Enter email addresses (one per line) in the text area
3. Click **Add Members**

> Note: Users must have Microreact accounts with matching email addresses.

#### Removing Members

1. Enter email addresses to remove (one per line)
2. Click **Remove Members**

#### Listing Members

1. Click **List Members** in the Team Management section
2. Current members will display in the output area

### Automatic Team Sharing

When uploading a project, you can automatically share it with your team:

1. Configure a Team ID in **Settings**
2. Go to the **Microreact** tab
3. Check **Share with team after upload**
4. Select the **Team role**:
   - **Viewer**: Can view but not edit
   - **Editor**: Can make changes to the project
   - **Manager**: Full control including sharing permissions
5. Upload your project

---

## Viewing Projects

### From the Application

After uploading, your project URL appears in the **Microreact Visualization** panel:

- Click the URL link to open in a new browser tab
- Click **Load in Viewer Below** to view embedded in the application

### From Microreact.org

1. Go to [microreact.org](https://microreact.org)
2. Sign in to your account
3. Click **My Projects** to see all your uploaded projects
4. Click any project to open it

### Project History

The **Previous Microreact Projects** table shows all projects uploaded during the current session, with clickable links to each.

---

## Troubleshooting

### "Please configure your Microreact API token"

- Go to **Settings** tab and enter your API token
- Make sure to click **Save Token** after pasting

### "Token format may be incorrect"

- Microreact tokens are JWT format and start with `eyJ`
- Make sure you copied the complete token from your account settings
- Check for extra spaces or line breaks

### "Upload failed" errors

Common causes:
- **Invalid token**: Token may have expired or been regenerated
- **Network issues**: Check your internet connection
- **Invalid project data**: Ensure your tree and metadata are valid

### Team operations fail

- Verify your API token has not expired
- Ensure the Team ID is correct (case-sensitive)
- Check that email addresses are valid Microreact accounts

### Project doesn't show tree/map

- Ensure your metadata has an `id` column matching tree tip labels
- For maps, the app now uses **ISO 3166-1 alpha-2 country codes** (e.g., "NO" for Norway, "SE" for Sweden)
- The app automatically converts country names to ISO codes for mapping
- If you have latitude/longitude columns, those will be used instead
- Check that the tree file is valid Newick format

### Countries not showing on map

The app includes a mapping of common country names to ISO codes. If your country isn't appearing:
- Check the country name in your data matches a known name (e.g., "Norway", "Sweden", "Germany")
- Some variations are supported: "Norge" → "NO", "Sverige" → "SE", "Tyskland" → "DE"
- Unknown countries will not appear on the map

---

## API Reference

The following R functions are available in `R/microreact_api.R` for programmatic use:

### Project Management

```r
# Upload a project from JSON
result <- microreact_create_project(token, project_json)
# Returns: list(success, data$id, data$url, error)

# Upload from a .microreact file
result <- microreact_create_project_from_file(token, file_path)

# Update an existing project
result <- microreact_update_project(token, project_id, project_json)
```

### Team Management

```r
# Create a new team
result <- microreact_create_team(token, "Team Name")
# Returns: list(success, data$id, error)

# Add members to a team
result <- microreact_add_team_members(token, team_id, c("user@example.com"))

# List team members
result <- microreact_list_team_members(token, team_id)

# Remove members from a team
result <- microreact_remove_team_members(token, team_id, c("user@example.com"))
```

### Project Sharing

```r
# Share a project with a team
result <- microreact_share_with_team(token, team_id, project_id, role = "viewer")
# role can be: "viewer", "editor", or "manager"

# Remove team access to a project
result <- microreact_unshare_from_team(token, team_id, project_id)
```

### Utility Functions

```r
# Validate token format
is_valid <- microreact_validate_token_format(token)

# Build project JSON from metadata and tree
json <- build_microreact_json(metadata_csv, tree_newick, project_name, description)
```

---

## External Resources

- [Microreact Documentation](https://docs.microreact.org/)
- [Microreact API - Creating Projects](https://docs.microreact.org/api/creating-projects)
- [Microreact API - Access Tokens](https://docs.microreact.org/api/access-tokens)
- [Microreact GitHub](https://github.com/microreact)

---

*Last updated: January 2026*
